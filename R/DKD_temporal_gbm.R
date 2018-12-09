setwd("~/dkd/DKD_PM")

rm(list=ls()); gc()
source("../helper_functions.R")
source("./utils.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"xgboost"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
))

# Load data
load("./data/DKD_heron_facts_prep.Rdata")
load("./data/pat_episode.Rdata")

##global values
#hyper-parameter grid for xgboost
eval_metric<-"auc"
objective<-"binary:logistic"
grid_params_tree<-expand.grid(
                              max_depth=10,
                              # max_depth=10,
                              eta=c(0.02,0.01),
                              # eta=0.01,
                              min_child_weight=1,
                              subsample=0.8,
                              colsample_bytree=0.8, 
                              gamma=1
                              )

#==== evaluations
# time_iterv<-"3mth"
time_iterv<-"6mth"
# time_iterv<-"1yr"

# ep_unit<-90
ep_unit<-182.5
# ep_unit<-365.25

#different temporal data handling methods
type_seq<-c("non-temporal",
            "stack-temporal",
            "discrete-survival")
nfold<-5 #nfold for cv

#result holder
result_lst<-list()
bm_lst<-list()

for(k in seq_along(type_seq)){
  bm<-c()
  bm_nm<-c()

  type<-type_seq[k]
  
  start_k<-Sys.time()
  cat("...temporal data integrated by:",type,'\n')
  
  #==build training set
  start_k_i<-Sys.time()
  
  train_pt<-pat_episode %>% 
    filter(!(part_idx %in% c("H")) & ((DM_YR+(episode/(365.25/ep_unit))) < 2016))

  if(type=="non-temporal"){
    train<-get_most_recent(train_pt,fact_stack,ep_unit)
    train_pt %<>%
      semi_join(train$y_train,by="PATIENT_NUM")

  }else if(type=="stack-temporal"){
    train<-get_stack_temporal(train_pt,fact_stack,ep_unit)
    train_pt %<>%
      semi_join(train$y_train,by=c("PATIENT_NUM","episode"))

  }else if(type=="discrete-survival"){
    train<-get_surv_temporal(train_pt,fact_stack,ep_unit)
    train_pt %<>%
      semi_join(train$y_train,by=c("PATIENT_NUM","episode"))

  }else{
    stop("temporal data processing method hasn't been implemented!")
  }
  dtrain<-xgb.DMatrix(data=train$x_train,
                      label=train$y_train$DKD_IND_episode)
  
  lapse_k_i<-Sys.time()-start_k_i
  cat("......finish handling temporal data in",lapse_k_i,units(lapse_k_i),"\n")
  bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
  bm_nm<-c(bm_nm,"build temporal traning set")
  
  #==get indices for k folds
  folds<-list()
  if(type %in% c("non-temporal")) {
    train_dt<-train_pt %>%
      dplyr::select(PATIENT_NUM,part_idx) %>% unique %>%
      dplyr::mutate(row_idx = 1:n())
  }else{
    train_dt<-train_pt %>%
      dplyr::mutate(row_idx = 1:n())
  }
  
  for(fd in seq_len(nfold)){
    fd_df<-train_dt %>% filter(part_idx==as.character(fd)) %>%
      dplyr::select(row_idx)
    folds[[fd]]<-fd_df$row_idx
  }

  #==tune hyperparameter
  start_k_i<-Sys.time()
  
  verb<-TRUE
  bst_grid<-c()
  bst_grid_cv<-c()
  metric_name<-paste0("test_", eval_metric,"_mean")
  metric_sd_name<-paste0("test_", eval_metric,"_std")
  # grid_params<-grid_params_reg
  grid_params<-grid_params_tree
  
  for(i in seq_len(dim(grid_params)[1])){
    start_i<-Sys.time()
    param<-as.list(grid_params[i,])
    param$scale_pos_weight=mean(train$y_train$DKD_IND_episode) #inbalance sampling
    # param$scale_pos_weight=1 #balance sampling
    
    if(type=="discrete-survival"){
      dtrain_b<-xgb.DMatrix(data=train$x_train_base,
                            label=train$y_train$DKD_IND_episode)
      
      #force an initial tree with only period info
      baseline<-xgb.cv(param,
                       dtrain_b,
                       objective = objective,
                       metrics = eval_metric,
                       maximize = TRUE,
                       nrounds=1,
                       # nfold = 5,
                       folds = folds,
                       prediction = T)
      
      setinfo(dtrain,"base_margin",baseline$pred)
    }
    
    bst <- xgb.cv(param,
                  dtrain,
                  objective = objective,
                  metrics = eval_metric,
                  maximize = TRUE,
                  nrounds=2000,
                  # nfold = 5,
                  folds = folds,
                  early_stopping_rounds = 100,
                  print_every_n = 100,
                  prediction = T) #keep cv results
    
    bst_grid<-rbind(bst_grid, cbind(grid_params[i,],
                                    metric=max(bst$evaluation_log[[metric_name]]),
                                    steps=which(bst$evaluation_log[[metric_name]]==max(bst$evaluation_log[[metric_name]]))[1]))
    
    bst_grid_cv<-cbind(bst_grid_cv,bst$pred)
    
    if(verb){
      cat('finished train case:',paste0(paste0(c(colnames(grid_params),"scale_pos_weight"),"="),param,collapse="; "),
          'in',Sys.time()-start_i,units(Sys.time()-start_i),"\n")
      start_i<-Sys.time()
    }
  }
  hyper_param<-bst_grid[which.max(bst_grid$metric),]
  valid_cv<-bst_grid_cv[,which.max(bst_grid$metric)]
  
  valid_cv<-data.frame(PATIENT_NUM = row.names(train$x_train),
                       episode = train$y_train$episode,
                       valid_type = 'cv',
                       pred = valid_cv,
                       real = getinfo(dtrain,"label"),
                       method = type,
                       stringsAsFactors = F)
  
  lapse_k_i<-Sys.time()-start_k_i
  cat("......finish tunning hyperparameter in",lapse_k_i,units(lapse_k_i),"\n")
  bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
  bm_nm<-c(bm_nm,"tune hyperparameters")
  
  #==validation
  start_k_i<-Sys.time()
  
  col_tr<-train$col_tr %>% arrange(CONCEPT_CD)
  
  test_pt<-pat_episode %>% 
    filter(part_idx %in% c("H") | ((DM_YR+(episode/(365.25/ep_unit))) >= 2016))
  
  test<-get_test_temporal(fact_stack,test_pt,col_tr)
  
  dtest<-xgb.DMatrix(data=test$x_test,
                     label=test$y_test$DKD_IND_episode)
  
  if(type=="discrete-survival"){
    trb<-train$x_train_base
    tsb<-test$x_test_base
    inters<-intersect(colnames(trb),colnames(tsb))
    
    dtrain_b<-xgb.DMatrix(data=trb[,which(colnames(trb) %in% inters)],
                          label=train$y_train$DKD_IND_episode)

    baseline<-xgb.train(data=dtrain_b,
                        max_depth=hyper_param$max_depth,
                        maximize = TRUE,
                        eta=hyper_param$eta,
                        nrounds=1,
                        eval_metric="auc",
                        objective="binary:logistic")
    
    base_margin_t<-predict(baseline,newdata=tsb[,which(colnames(tsb) %in% inters)],outputmargin=TRUE)
    
    #update base margin by at different interval
    setinfo(dtest,"base_margin",base_margin_t)
  }

  watchlist<-list(test=dtest)
  xgb_tune<-xgb.train(data=dtrain,
                      max_depth=hyper_param$max_depth,
                      maximize = TRUE,
                      eta=hyper_param$eta,
                      nrounds=hyper_param$steps,
                      watchlist=watchlist,
                      eval_metric="auc",
                      objective="binary:logistic",
                      print_every_n = 100)
  
  valid<-data.frame(PATIENT_NUM = as.character(test$y_test$PATIENT_NUM),
                    episode = test$y_test$episode,
                    valid_type = test$y_test$part_idx,
                    pred = predict(xgb_tune,dtest),
                    real = test$y_test$DKD_IND_episode,
                    method = type,
                    stringsAsFactors = F)
  
  lapse_k_i<-Sys.time()-start_k_i
  cat("......finish validating in",lapse_k_i,units(lapse_k_i),"\n")
  bm<-c(bm,lapse_k_i)
  bm_nm<-c(bm_nm,"validation")
  
  lapse_k<-Sys.time()-start_k
  cat("...finish evaluating temporal integraion method:",type,"in",lapse_k,units(lapse_k),"\n")
  bm<-c(bm,lapse_k_i)
  bm_nm<-c(bm_nm,"end method loop")
  
  result_lst[[type]]<-valid %>% bind_rows(valid_cv)
  bm_lst[[type]]<-data.frame(step=bm_nm,time=bm)
}

save(result_lst,file=paste0("./data/",time_iterv,"_","temporal_gbm.Rdata"))
save(bm_lst,file=paste0("./data/",time_iterv,"_","temporal_gbm_bm.Rdata"))

rm(list=ls())
gc()


