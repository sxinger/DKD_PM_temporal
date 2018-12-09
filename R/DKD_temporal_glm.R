setwd("~/dkd/DKD_PM")

rm(list=ls()); gc()
source("../helper_functions.R")
source("./utils.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"h2o"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
))

# Load data
load("./data/DKD_heron_facts_prep.Rdata")
# load("./data/pat_episode.Rdata")
load("./data/pat_episode2.Rdata")

##global values
#hyper-parameter grid for xgboost
hyper_params<-list(alpha=c(1,0.8,0.5,0.3))

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

#initialize h2o
h2o.init(nthreads=-1)

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
    
    dkd_fold<-train$y_train %>%
      inner_join(pat_episode %>%
                   dplyr::select(PATIENT_NUM,fold) %>% unique,
                 by=c("PATIENT_NUM"))
    
  }else if(type=="stack-temporal"){
    train<-get_stack_temporal(train_pt,fact_stack,ep_unit)
    train_pt %<>%
      semi_join(train$y_train,by=c("PATIENT_NUM","episode"))
    
    dkd_fold<-train$y_train %>% 
      inner_join(pat_episode %>%
                   dplyr::select(PATIENT_NUM,episode,fold) %>% unique,
                 by=c("PATIENT_NUM","episode"))
    
  }else if(type=="discrete-survival"){
    train<-get_surv_temporal(train_pt,fact_stack,ep_unit)
    train_pt %<>%
      semi_join(train$y_train,by=c("PATIENT_NUM","episode"))
    
    dkd_fold<-train$y_train %>% 
      inner_join(pat_episode %>%
                   dplyr::select(PATIENT_NUM,episode,fold) %>% unique,
                 by=c("PATIENT_NUM","episode"))
    
  }else{
    stop("temporal data processing method hasn't been implemented!")
  }
  
  
  lapse_k_i<-Sys.time()-start_k_i
  cat("......finish handling temporal data in",lapse_k_i,units(lapse_k_i),"\n")
  bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
  bm_nm<-c(bm_nm,"build temporal traning set")
  
  #==tune hyperparameter
  start_k_i<-Sys.time()
  
  #combine predictors and target (as required by h2o)
  train_mt<-cbind(fold_assign=dkd_fold$fold,
                  DKD_IND=train$y_train$DKD_IND_episode,
                  train$x_train)
  
  predb_idx<-which(grepl('^(ep_)+',colnames(train_mt)))
  pred_idx<-which(!colnames(train_mt) %in% c("fold_assign","DKD_IND"))
  target_idx<-which(colnames(train_mt)=="DKD_IND")
  fold_idx<-which(colnames(train_mt)=="fold_assign")
  
  #h2o bug: special characters in column names will result in an additional row
  col_encode<-data.frame(col_name=colnames(train_mt),
                         col_code=paste0("C",1:ncol(train_mt)))
  colnames(train_mt)<-col_encode$col_code
  train_h2o<-as.h2o(as.matrix(train_mt)) #column order doesn't change
  #h2o truncate column space at 100 when apply as.h2o
  #h2o.splitFrame would force index and print out all columns

  if(type=="discrete-survival"){
    #force an initial model including only period info
    baseline<-h2o.glm(x=predb_idx,
                      y=target_idx,   
                      training_frame=train_h2o,
                      family="binomial",
                      solver="COORDINATE_DESCENT",   #same optimization method as glmnet
                      alpha=0,                       #ridge regression
                      fold_column = paste0("C",fold_idx),
                      lambda_search=TRUE,
                      early_stopping = TRUE,
                      # missing_values_handling="Skip",
                      remove_collinear_columns=TRUE,
                      keep_cross_validation_predictions =T
                     )
    
    log_os<-h2o.getFrame(baseline@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
    os<-log(as.data.frame(log_os)$p1/(1-as.data.frame(log_os)$p1))
    
    train_h2o<-h2o.cbind(train_h2o,as.h2o(data.frame(offset=os)))
    
    alpha_grid<-h2o.grid(x=pred_idx,
                         y=target_idx,  
                         training_frame=train_h2o,
                         algorithm = "glm",
                         grid_id =paste0("elastnet_alpha_grid",k),
                         family="binomial",
                         solver="COORDINATE_DESCENT",   #same optimization method as glmnet
                         fold_column = paste0("C",fold_idx),
                         offset_column = "offset",
                         lambda_search=TRUE,
                         early_stopping = TRUE,
                         # missing_values_handling="Skip",
                         remove_collinear_columns=TRUE,
                         keep_cross_validation_predictions =T,
                         hyper_params = hyper_params,
                         search_criteria = list(strategy = "Cartesian")
                         )
  }else{
    alpha_grid<-h2o.grid(x=pred_idx,
                         y=target_idx,  
                         training_frame=train_h2o,
                         algorithm = "glm",
                         grid_id = paste0("elastnet_alpha_grid",k),
                         family="binomial",
                         solver="COORDINATE_DESCENT",   #same optimization method as glmnet
                         fold_column = paste0("C",fold_idx),
                         lambda_search=TRUE,
                         early_stopping = TRUE,
                         # missing_values_handling="Skip",
                         remove_collinear_columns=TRUE,
                         keep_cross_validation_predictions =T,
                         hyper_params = hyper_params,
                         search_criteria = list(strategy = "Cartesian")
                         )
  }

  # alpha search for elastic net
  sortedGrid<-h2o.getGrid(paste0("elastnet_alpha_grid",k),sort_by = "auc",decreasing = TRUE)
  alpha_opt_model1<-h2o.getModel(sortedGrid@model_ids[[1]])
  alpha_opt<-alpha_opt_model1@parameters$alpha
  if(alpha_opt < 1){
    alpha_opt_model<-alpha_opt_model1
    lasso_model<-h2o.getModel(alpha_grid@model_ids[[length(hyper_params$alpha)]]) #models are queued
  }else{
    alpha_opt_model<-h2o.getModel(sortedGrid@model_ids[[2]])
    lasso_model<-alpha_opt_model1
    alpha_opt<-alpha_opt_model@parameters$alpha
  }
  
  lasso_cv<-h2o.getFrame(lasso_model@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
  elastnet_cv<-h2o.getFrame(alpha_opt_model@model[["cross_validation_holdout_predictions_frame_id"]][["name"]])
  
  valid_cv<-data.frame(PATIENT_NUM = row.names(train$x_train),
                       episode = train$y_train$episode,
                       valid_type = 'cv',
                       pred_lasso = as.data.frame(lasso_cv)$p1,
                       pred_elastnet = as.data.frame(elastnet_cv)$p1,
                       alpha=alpha_opt,
                       real = train$y_train$DKD_IND_episode,
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
  
  test<-get_test_temporal(fact_stack,pat_episode,col_tr)
  
  test_mt<-cbind(fold_assign=NA,
                 DKD_IND=test$y_test$DKD_IND_episode,
                 test$x_test)
  if(!all(colnames(test_mt)==col_encode$col_name)){
    stop("predictors of testing set not aligned with training set")
  }
  colnames(test_mt)<-col_encode$col_code
  test_h2o<-as.h2o(as.matrix(test_mt))
  
  if(type == "discrete-survival"){
    #combine predictors and target (as required by h2o)
    valid_os<-predict(baseline,test_h2o)
    logit_os<-as.data.frame(valid_os)$p1/(1-as.data.frame(valid_os)$p1)
    
    test_h2o<-h2o.cbind(test_h2o,as.h2o(data.frame(offset=logit_os)))
  }
  
  lasso_h<-predict(lasso_model,test_h2o)
  elastnet_h<-predict(alpha_opt_model,test_h2o)
    
  valid<-data.frame(PATIENT_NUM = as.character(test$y_test$PATIENT_NUM),
                    episode = test$y_test$episode,
                    valid_type = test$y_test$part_idx,
                    pred_lasso = as.data.frame(lasso_h)$p1,
                    pred_elastnet = as.data.frame(elastnet_h)$p1,
                    alpha = alpha_opt,
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

save(result_lst,file=paste0("./data/",time_iterv,"_","temporal_glm.Rdata"))
save(bm_lst,file=paste0("./data/",time_iterv,"_","temporal_glm_bm.Rdata"))

h2o.shutdown(prompt = FALSE)

rm(list=ls())
gc()

