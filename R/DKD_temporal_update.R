## landmark boosting ##
rm(list=ls()); gc()

# setwd("~/proj_dkd/DKD_PM_wip")

source("./R/util.R")
source("./R/newXY.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"xgboost"
                     ,"h2o"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
                  ))

# Load data
X_long<-readRDS("./data2/X_long.rda") %>%
  anti_join(readRDS("./data2/pat_T1DM.rda"),by="PATIENT_NUM")

pat_tbl<-readRDS("./data2/pat_episode2.rda") %>%
  anti_join(readRDS("./data2/pat_T1DM.rda"),by="PATIENT_NUM")

# time_iterv<-"3mth"
# time_iterv<-"6mth"
time_iterv<-"1yr"

# ep_unit<-90
# ep_unit<-182.5
ep_unit<-365.25

#period
period<-(365.25/ep_unit)*5
# period<-(365/ep_unit)*8

pat_episode<-pat_tbl %>%
  dplyr::mutate(episode = floor(as.numeric(DAY_SINCE)/ep_unit))

X_long %<>% dplyr::filter(episode_x < period)

# number of fold for cv
nfold<-5 #nfold for cv

# temporal handling method
type<-"dynamic-temporal"


model_yr<-list()
model_roc<-list()
bm_yr<-list()

for(yr in 0:4){
  for(past in c(T,F)){
    if(past){
      partition<-pat_episode %>%
        semi_join(X_long, by = "PATIENT_NUM") %>%
        group_by(PATIENT_NUM, episode) %>%
        top_n(n=1,wt=DKD_IND_additive) %>% ungroup %>%
        dplyr::select(PATIENT_NUM,episode,part_idx,DKD_IND_additive) %>%
        dplyr::mutate(PATIENT_NUM2=PATIENT_NUM,
                      episode2=episode) %>%
        unite("PATIENT_NUM_ep",c("PATIENT_NUM2","episode2")) %>%
        arrange(PATIENT_NUM_ep) %>% unique 
    }else{
      partition<-pat_episode %>%
        semi_join(X_long, by = "PATIENT_NUM") %>%
        group_by(PATIENT_NUM, episode) %>%
        top_n(n=1,wt=DKD_IND_additive) %>% ungroup %>%
        dplyr::select(PATIENT_NUM,episode,part_idx,DKD_IND_additive) %>%
        dplyr::mutate(PATIENT_NUM2=PATIENT_NUM,
                      episode2=episode) %>%
        unite("PATIENT_NUM_ep",c("PATIENT_NUM2","episode2")) %>%
        arrange(PATIENT_NUM_ep) %>% unique %>%
        filter(!(part_idx %in% c(5)))
    }
    
    start_yr<-Sys.time()
    bm<-c()
    bm_nm<-c()
    
    #==build training set
    start_i<-Sys.time()
    
    train_pt<-partition %>% filter(episode<=yr)
    
    Xy_tr<-get_dsurv_temporal(train_pt,X_long) #
    
    dtrain<-xgb.DMatrix(data=Xy_tr$X_ep,
                        label=Xy_tr$y_ep$DKD_IND_additive)
    
    if(yr>0){
      os<-data.frame(PATIENT_NUM_ep=row.names(Xy_tr$X_ep),stringsAsFactors = F) %>%
        left_join(model_roc[[paste0((yr-1),"yr_since_dm_TRUE")]] %>% 
                    dplyr::select(PATIENT_NUM,pred),
                  by=c("PATIENT_NUM_ep"="PATIENT_NUM")) %>%
        dplyr::mutate(pred_imp=mean(pred,na.rm=T)) %>%
        dplyr::mutate(offset=ifelse(is.na(pred),
                                    log(pred_imp/(1-pred_imp)),
                                    log(pred/(1-pred))))
      
      setinfo(dtrain,"base_margin",os$offset)
    }
    
    lapse_i<-Sys.time()-start_i
    cat("...finish handling temporal data in",lapse_i,units(lapse_i),"\n")
    bm<-c(bm,paste0(round(lapse_i,4),units(lapse_i)))
    bm_nm<-c(bm_nm,"build temporal traning set")
    
    #==get indices for k folds
    folds<-list()
    train_dt<-train_pt %>%
      dplyr::select(PATIENT_NUM,part_idx) %>% unique %>%
      dplyr::mutate(row_idx = 1:n())
    
    for(fd in seq_len(max(train_dt$part_idx))){
      fd_df<-train_dt %>% filter(part_idx==as.character(fd)) %>%
        dplyr::select(row_idx)
      folds[[fd]]<-fd_df$row_idx
    }
    
    #==tune hyperparameter
    #hyper-parameter grid for xgboost
    eval_metric<-"auc"
    objective<-"binary:logistic"
    grid_params_tree<-expand.grid(
      max_depth=8,
      # max_depth=c(2,8),
      eta=0.02,
      # eta=c(0.02,0.01),
      min_child_weight=1,
      subsample=0.8,
      colsample_bytree=0.8, 
      gamma=1
    )
    
    start_i<-Sys.time()
    
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
      # param$scale_pos_weight=mean(train$y_train$DKD_IND_additive) #inbalance sampling
      param$scale_pos_weight=1 #balance sampling
      
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
    
    valid_cv<-data.frame(PATIENT_NUM = row.names(Xy_tr$X_ep),
                         pred = bst_grid_cv[,which.max(bst_grid$metric)],
                         real = getinfo(dtrain,"label"),
                         stringsAsFactors = F)
    
    lapse_i<-Sys.time()-start_i
    cat("......finish tunning hyperparameter in",lapse_i,units(lapse_i),"\n")
    bm<-c(bm,paste0(round(lapse_i,4),units(lapse_i)))
    bm_nm<-c(bm_nm,"tune hyperparameters")
    
    ##==re-train model on all
    start_i<-Sys.time()
    
    xgb_tune<-xgb.train(data=dtrain,
                        max_depth=hyper_param$max_depth,
                        maximize = TRUE,
                        eta=hyper_param$eta,
                        nrounds=hyper_param$steps,
                        # watchlist=watchlist,
                        eval_metric="auc",
                        objective="binary:logistic",
                        print_every_n = 100)
    
    lapse_i<-Sys.time()-start_i
    cat("......finish validating in",lapse_i,units(lapse_i),"\n")
    bm<-c(bm,paste0(round(lapse_i,4),units(lapse_i)))
    bm_nm<-c(bm_nm,"validation")
    
    ##==save results
    model_yr[[paste0(yr,"yr_since_dm_",past)]]<-xgb_tune
    model_roc[[paste0(yr,"yr_since_dm_",past)]]<-valid_cv
    bm_yr[[paste0(yr,"yr_since_dm",past)]]<-data.frame(step=bm_nm,time=bm,stringsAsFactors = F)
  }
}

out<-list(model_yr = model_yr,
          model_roc = model_roc,
          bm_yr = bm_yr)

save(out,file=paste0("./data2/",time_iterv,"_",type,"_gbm_model4.Rdata"))
