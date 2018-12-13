## validations ##

setwd("~/proj_dkd/DKD_PM_temporal")

rm(list=ls()); gc()
source("../helper_functions.R")
source("./newXY.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"xgboost"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
                  ))

##=================collect validation results=================
validation<-list()

# Load patient data
load("./data2/pat_episode2.Rdata")

ep_unit<-365.25
pat_episode<-pat_tbl %>%
  dplyr::mutate(episode = floor(as.numeric(DAY_SINCE)/ep_unit))

partition_ts<-pat_episode %>%
  dplyr::select(PATIENT_NUM,part_idx,DAY_SINCE,episode,DKD_IND_additive) %>%
  group_by(PATIENT_NUM,part_idx) %>%
  top_n(n=1,wt=DAY_SINCE) %>% ungroup %>%
  unique %>% filter(part_idx %in% c(5))

# temporal handling types
type_seq<-c(
  "non-temporal"
 ,"stack-temporal"
 ,"discrt-surv-temporal"
 ,"grp-temporal"
 ,"dynamic-temporal"
)

# model<-"glm"
model<-"gbm"

# time_iterv<-"3mth"
# time_iterv<-"6mth"
time_iterv<-"1yr"

####dynamic-temporal####
i<-5
load(paste0("./data2/",time_iterv,"_",type_seq[i],"_gbm_model2.Rdata"))
load("./data/X_long.Rdata")

#subset testing set
X_long %<>%
  semi_join(partition_ts %>% dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)),
            by="PATIENT_NUM")

valid_stack5<-c()
for(j in 0:4){
  start<-Sys.time()
  
  #get predictors for test set
  test<-X_long %>% filter(episode_x < j) %>%
    dplyr::mutate(episode=j) %>%
    unite("PATIENT_NUM_ep",c("PATIENT_NUM","episode")) %>%
    unite("VARIABLE_ep",c("VARIABLE","episode_x")) %>%
    dplyr::select(PATIENT_NUM_ep, VARIABLE_ep, NVAL_NUM) %>%
    unique 
  
  test1<-test %>%
    semi_join(data.frame(VARIABLE_ep = out$model_yr[[paste0(j,"yr_since_dm_FALSE")]]$feature_names,
                         stringsAsFactors = F),
              by="VARIABLE_ep")
  
  test2<-test %>%
    semi_join(data.frame(VARIABLE_ep = out$model_yr[[paste0(j,"yr_since_dm_TRUE")]]$feature_names,
                         stringsAsFactors = F),
              by="VARIABLE_ep")
    
  
  #collect additional features required for training model
  x_add1<-data.frame(VARIABLE_ep = out$model_yr[[paste0(j,"yr_since_dm_FALSE")]]$feature_names,
                    stringsAsFactors = F) %>%
    anti_join(data.frame(VARIABLE_ep = as.character(unique(test1$VARIABLE_ep)),
                         stringsAsFactors = F),
              by="VARIABLE_ep")
  
  x_add2<-data.frame(VARIABLE_ep = out$model_yr[[paste0(j,"yr_since_dm_TRUE")]]$feature_names,
                     stringsAsFactors = F) %>%
    anti_join(data.frame(VARIABLE_ep = as.character(unique(test2$VARIABLE_ep)),
                         stringsAsFactors = F),
              by="VARIABLE_ep")
  
  #introduce psuedo-obs to align with training model
  if(nrow(x_add1)>0){
    test1 %<>%
      bind_rows(data.frame(PATIENT_NUM_ep = rep("0_0",nrow(x_add1)),
                           VARIABLE_ep = x_add1$VARIABLE_ep,
                           NVAL_NUM = 0,
                           stringsAsFactors = F))
  }
  
  if(nrow(x_add2)>0){
    test2 %<>%
      bind_rows(data.frame(PATIENT_NUM_ep = rep("0_0",nrow(x_add2)),
                           VARIABLE_ep = x_add2$VARIABLE_ep,
                           NVAL_NUM = 0,
                           stringsAsFactors = F))
  }
  
  test1 %<>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM_ep",
                          variable="VARIABLE_ep",
                          val="NVAL_NUM")
  
  test2 %<>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM_ep",
                          variable="VARIABLE_ep",
                          val="NVAL_NUM")
  
  if(nrow(x_add1)>0){
    test1<-test1[-1,]
  }
  
  if(nrow(x_add2)>0){
    test2<-test2[-1,]
  }
  
  #check alignment
  if(!all(colnames(test1)==colnames(out$model_yr[[paste0(j,"yr_since_dm_FALSE")]]$feature_names))){
    warning("predictors don't align!")
  }
  if(!all(colnames(test2)==colnames(out$model_yr[[paste0(j,"yr_since_dm_TRUE")]]$feature_names))){
    warning("predictors don't align!")
  }

  #make prediction on validation set
  wt<-pROC::auc(out$model_roc[[paste0(j,"yr_since_dm_TRUE")]]$real,
                out$model_roc[[paste0(j,"yr_since_dm_TRUE")]]$pred)
  
  valid<-data.frame(PATIENT_NUM_ep = row.names(test1),
                    episode = j,
                    pred1 = predict(out$model_yr[[paste0(j,"yr_since_dm_FALSE")]],
                                   newdata=test1),
                    pred2 = predict(out$model_yr[[paste0(j,"yr_since_dm_TRUE")]],
                                    newdata=test2),
                    wt=as.numeric(wt),
                    stringsAsFactors = F) %>%
    dplyr::mutate(pred=pred1)
  
  if(j>0){
    valid %<>%
      dplyr::mutate(PATIENT_NUM_ep_last = PATIENT_NUM_ep) %>%
      separate("PATIENT_NUM_ep_last",c("PATIENT_NUM","episode")) %>%
      dplyr::mutate(episode=as.numeric(episode)) %>%
      dplyr::mutate(episode_last = episode-1) %>%
      unite("PATIENT_NUM_ep_last",c("PATIENT_NUM","episode_last")) %>%
      left_join(valid_stack5 %>% filter(episode==j-1) %>%
                  dplyr::select(PATIENT_NUM_ep,pred2) %>%
                  dplyr::rename(pred_last=pred2),
                by=c("PATIENT_NUM_ep_last"="PATIENT_NUM_ep")) %>%
      dplyr::mutate(pred_imp=mean(pred_last,na.rm=T)) %>%
      dplyr::mutate(offset=ifelse(is.na(pred_last),
                                  log(pred_imp/(1-pred_imp)),
                                  log(pred_last/(1-pred_last)))) %>%
      dplyr::mutate(pred_linear=log(pred1/(1-pred1))+offset) %>%
      dplyr::mutate(pred = exp(pred_linear)/(1+exp(pred_linear))) %>%
      dplyr::select(PATIENT_NUM_ep,episode,pred1,pred2,wt,pred)
  }
  
  valid_stack5 %<>% bind_rows(valid)
  
  lapse<-Sys.time()-start
  cat("finish validate episode",j,"in",lapse,units(lapse),".\n")
}

partition_ts %<>%
  dplyr::mutate(PATIENT_NUM = as.character(PATIENT_NUM),
                real=DKD_IND_additive) %>%
  dplyr::select(PATIENT_NUM,episode,real,part_idx)


valid_stack5 %<>%
  dplyr::mutate(PATIENT_NUM = gsub("_.*","",PATIENT_NUM_ep)) %>%
  left_join(partition_ts,by=c("PATIENT_NUM","episode")) %>%
  dplyr::filter(!((part_idx=="H") && (episode > 1))) %>%
  dplyr::mutate(real_ref=real) %>%
  group_by(PATIENT_NUM) %>%
  tidyr::fill(.,c("real_ref","part_idx"),.direction="up") %>%
  ungroup %>% dplyr::filter(!is.na(real_ref)) %>%
  replace_na(list(real=0)) %>%
  group_by(PATIENT_NUM,part_idx) %>%
  arrange(episode) %>%
  dplyr::mutate(pred2=cumsum(pred),
                pred3=cummax(pred),
                pred4=cumsum(1/pred),
                wt2=cumsum(wt),
                cum_row=1:n()) %>%
  ungroup %>%
  dplyr::mutate(pred4=cum_row/pred4,
                pred_wt=pred2*wt/wt2,
                pred_wt2=pred2/cum_row) %>%
  group_by(episode) %>%
  dplyr::select(PATIENT_NUM,episode,pred,pred2,pred3,pred4,
                pred_wt,pred_wt2,
                real,part_idx) %>%
  dplyr::mutate(temporal=type_seq[i]) %>%
  ungroup %>%
  semi_join(pat_episode %>% dplyr::select(PATIENT_NUM,episode) %>% 
              unique %>% dplyr::mutate(PATIENT_NUM=as.character(PATIENT_NUM)),
            by=c("PATIENT_NUM","episode"))

load("./data2/valid_stack5.Rdata")
valid_stack5 %>%
  semi_join(pat_episode %>% dplyr::select(PATIENT_NUM,episode) %>% 
              unique %>% dplyr::mutate(PATIENT_NUM = as.character(PATIENT_NUM)),
            by=c("PATIENT_NUM","episode")) %>%
  dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)) %>%
  group_by(episode,part_idx,temporal) %>%
  dplyr::summarise(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM*real))-1,
                   dkd_rate = round((length(unique(PATIENT_NUM*real))-1)/length(unique(PATIENT_NUM)),2),
                   low95=pROC::ci.auc(real,pred)[1],
                   auc=pROC::ci.auc(real,pred)[2],
                   up95=pROC::ci.auc(real,pred)[3]
                   ) %>%
  ungroup %>% arrange(part_idx, episode) %>%
  View

validation[[paste0("valid_stack",i)]]<-valid_stack5 %>% ungroup

####grp-temporal####
i<-4
load(paste0("./data/",time_iterv,"_",type_seq[i],"_gbm_model.Rdata"))
load("./data/X_long.Rdata")

valid_stack4<-c()
for(j in c(0,seq_len(max(partition_ts$episode)))){
  start_j<-Sys.time()
  
  #get predictors for test set
  test<-X_long %>% 
    left_join(partition_ts %>%
                dplyr::select(PATIENT_NUM,episode) %>% 
                dplyr::filter(episode == j) %>% unique,
              by="PATIENT_NUM") %>%
    dplyr::filter(episode > episode_x) %>%
    dplyr::select(PATIENT_NUM, VARIABLE, episode_x, NVAL_NUM) %>%
    unique
  
  #standardization
  test %<>%
    inner_join(out$x_prep$Tr_std,by="VARIABLE") %>%
    dplyr::mutate(NVAL_std=(NVAL_NUM-NVAL_mean)/NVAL_sd)
  
  #orthogonalization
  var_lst<-names(out$x_prep$Tr_ort)
  test_std<-c()
  
  for(var in var_lst){
    start<-Sys.time()
    
    #collect episodes used for training
    ep<-X_long %>% anti_join(partition_ts,by="PATIENT_NUM") %>%
      dplyr::filter(VARIABLE==var) %>% 
      dplyr::select(episode_x) %>% 
      arrange(episode_x) %>% unique %>%
      left_join(test %>% 
                  dplyr::filter(VARIABLE==var) %>%
                  dplyr::select(episode_x) %>%
                  arrange(episode_x) %>% unique %>%
                  dplyr::mutate(col_sel=T),
                by="episode_x") %>%
      replace_na(list(col_sel=F))
      
    
    #transform var
    if(all(ep$col_sel==F) | is.null(out$x_prep$Tr_ort[[var]])){
      test_std %>%
        bind_rows(test %>% dplyr::select(PATIENT_NUM, VARIABLE, NVAL_std) %>%
                    dplyr::filter(VARIABLE==var) %>%  unique)
    }else{
      var_trans<-test %>%
        dplyr::filter(VARIABLE==var) %>%
        dplyr::select(PATIENT_NUM,VARIABLE,episode_x,NVAL_std) %>%
        semi_join(ep %>% dplyr::filter(col_sel),by="episode_x") %>% #episode in test may not necessarily been selected by training
        unite("VARIABLE_ep",c("VARIABLE","episode_x")) %>%
        long_to_sparse_matrix(.,
                              id="PATIENT_NUM",
                              variable="VARIABLE_ep",
                              val="NVAL_std")
      var_trans<-var_trans %*% out$x_prep$Tr_ort[[var]][ep$col_sel,]
      colnames(var_trans)<-paste0(var,"_Comp",1:ncol(var_trans))
      
      #stack standardized features
      test_std %<>%
        bind_rows(data.frame(PATIENT_NUM = as.numeric(rep(row.names(var_trans),ncol(var_trans))),
                             VARIABLE    = rep(colnames(var_trans),each=nrow(var_trans)),
                             NVAL_std    = var_trans@x,
                             stringsAsFactors = F))
    }
    
    lapse<-Sys.time()-start
    cat("...",var,"transformed in",lapse,units(lapse),".\n")
  }
  
  #collect additional features required for training model
  x_add<-data.frame(VARIABLE = out$model$feature_names,
                    stringsAsFactors = F) %>%
    anti_join(data.frame(VARIABLE = as.character(unique(test_std$VARIABLE)),
                         stringsAsFactors = F),
              by="VARIABLE")
  
  #introduce psuedo-obs to align with training model
  test_std %<>%
    bind_rows(data.frame(PATIENT_NUM = rep(0,nrow(x_add)),
                         VARIABLE = x_add$VARIABLE,
                         NVAL_NUM = 0,
                         stringsAsFactors = F))
  
  test_std %<>%
    arrange(PATIENT_NUM) %>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM",
                          variable="VARIABLE",
                          val="NVAL_NUM")
  
  test_std<-test_std[-1,]
  
  #check alignment
  if(!all(colnames(test_std)==colnames(out$x_prep$X_ep))){
    warning("predictors don't align!")
  }
  
  #make prediction on validation set
  valid<-data.frame(PATIENT_NUM = as.numeric(row.names(test_std)),
                    pred = predict(out$model,newdata=test_std),
                    stringsAsFactors = F) %>%
    inner_join(partition_ts %>%
                 filter(episode == j) %>%
                 group_by(PATIENT_NUM,episode, part_idx) %>% 
                 top_n(n=1,wt=DKD_IND_additive) %>%
                 ungroup %>% dplyr::rename(real = DKD_IND_additive),
               by="PATIENT_NUM") %>%
    dplyr::select(PATIENT_NUM, episode, pred, real, part_idx) %>%
    dplyr::mutate(temporal = type_seq[i])
  
  valid_stack4 %<>% bind_rows(valid)
  
  lapse_j<-Sys.time()-start_j
  cat("finish predict at",j, "in",lapse_j,units(lapse_j),".\n")
}

#quick AUC check
valid_stack4 %>%
  group_by(episode,part_idx,temporal) %>%
  dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM*real))-1,
                   dkd_rate = round((length(unique(PATIENT_NUM*real))-1)/length(unique(PATIENT_NUM)),2),
                   low95=pROC::ci.auc(real,pred)[1],
                   auc=pROC::ci.auc(real,pred)[2],
                   up95=pROC::ci.auc(real,pred)[3]) %>%
  ungroup %>% arrange(part_idx, episode) %>%
  View

  
####discrt-surv-temporal####
i<-3
load(paste0("./data2/",time_iterv,"_",type_seq[i],"_gbm_model2.Rdata"))
load("./data/X_long.Rdata")

#subset testing set
X_long %<>% semi_join(partition_ts, by="PATIENT_NUM")

valid_stack3<-c()
for(j in 0:4){
  #get predictors for test set
  test<-X_long %>% filter(episode_x < j) %>%
    dplyr::mutate(episode=j) %>%
    unite("PATIENT_NUM_ep",c("PATIENT_NUM","episode")) %>%
    unite("VARIABLE_ep",c("VARIABLE","episode_x")) %>%
    dplyr::select(PATIENT_NUM_ep, VARIABLE_ep, NVAL_NUM) %>%
    unique
  
  test %<>%
    semi_join(data.frame(VARIABLE_ep = out$model$feature_names,
                         stringsAsFactors = F),
              by="VARIABLE_ep")
  
  #collect additional features required for training model
  x_add<-data.frame(VARIABLE_ep = out$model$feature_names,
                    stringsAsFactors = F) %>%
    anti_join(data.frame(VARIABLE_ep = as.character(unique(test$VARIABLE_ep)),
                         stringsAsFactors = F),
              by="VARIABLE_ep")
  
  #introduce psuedo-obs to align with training model
  if(nrow(x_add) > 0){
    test %<>%
      bind_rows(data.frame(PATIENT_NUM_ep = rep("0",nrow(x_add)),
                           VARIABLE_ep = x_add$VARIABLE_ep,
                           NVAL_NUM = 0,
                           stringsAsFactors = F))
  }

  test %<>%
    arrange(PATIENT_NUM_ep) %>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM_ep",
                          variable="VARIABLE_ep",
                          val="NVAL_NUM")
  if(nrow(x_add) > 0){
    test<-test[-1,] 
  }
  
  #check alignment
  if(!all(colnames(test)==colnames(out$x_prep$X_ep))){
    warning("predictors don't align!")
  }
  
  #make prediction on validation set
  valid<-data.frame(PATIENT_NUM = row.names(test),
                    episode=j,
                    pred = predict(out$model,newdata=test),
                    stringsAsFactors = F)
  
  valid_stack3 %<>% bind_rows(valid)
}


partition_ts %<>%
  dplyr::mutate(PATIENT_NUM = as.character(PATIENT_NUM),
                real=DKD_IND_additive) %>%
  dplyr::select(PATIENT_NUM,episode,real,part_idx)


valid_stack3 %<>%
  dplyr::mutate(PATIENT_NUM = gsub("_.*","",PATIENT_NUM)) %>%
  left_join(partition_ts,by=c("PATIENT_NUM","episode")) %>%
  dplyr::filter(!((part_idx=="H") && (episode > 1))) %>%
  dplyr::mutate(real_ref=real) %>%
  group_by(PATIENT_NUM) %>%
  tidyr::fill(.,c("real_ref","part_idx"),.direction="up") %>%
  ungroup %>% dplyr::filter(!is.na(real_ref)) %>%
  replace_na(list(real=0)) %>%
  group_by(PATIENT_NUM,part_idx) %>%
  arrange(episode) %>%
  dplyr::mutate(pred2=cumsum(pred),
                pred3=cummax(pred),
                pred4=cumsum(1/pred),
                cum_row=1:n()) %>%
  ungroup %>%
  dplyr::mutate(pred4=cum_row/pred4,
                pred_wt2=pred2/cum_row) %>%
  group_by(episode) %>%
  dplyr::mutate(pred5=scales::rescale(pred2,c(0.1,0.99))) %>%
  dplyr::select(PATIENT_NUM,episode,pred,pred2,pred3,pred4,pred_wt2,
                real,part_idx) %>%
  dplyr::mutate(temporal=type_seq[i]) %>%
  semi_join(pat_episode %>% dplyr::select(PATIENT_NUM,episode) %>% 
              unique %>% dplyr::mutate(PATIENT_NUM=as.character(PATIENT_NUM)),
            by=c("PATIENT_NUM","episode"))

#quick AUC check
valid_stack3 %>%
  dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)) %>%
  group_by(episode,part_idx,temporal) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM*real))-1,
                   dkd_rate = round((length(unique(PATIENT_NUM*real))-1)/length(unique(PATIENT_NUM)),2),
                   low95=pROC::ci.auc(real,pred)[1],
                   auc=pROC::ci.auc(real,pred)[2],
                   up95=pROC::ci.auc(real,pred)[3]) %>%
  ungroup %>% arrange(part_idx, episode) %>%
  View

validation[[paste0("valid_stack",i)]]<-valid_stack3 %>% ungroup

#####stack-temporal#####
i<-2
load(paste0("./data2/",time_iterv,"_",type_seq[i],"_gbm_model3.Rdata"))
load("./data/X_long.Rdata")

#subset testing set
X_long %<>%
  semi_join(partition_ts %>% dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)),
            by="PATIENT_NUM")

valid_stack2<-c()
for(j in 0:4){
  start<-Sys.time()
  
  model<-out$model_yr[[paste0(j,"yr_since_dm")]]
  feat_tr<-colnames(out$model_tr[[paste0(j,"yr_since_dm")]])
  
  #get predictors for test set
  test<-X_long %>% filter(episode_x < j) %>%
    dplyr::mutate(VARIABLE=paste0(VARIABLE,"_",episode_x)) %>%
    dplyr::select(PATIENT_NUM, VARIABLE, NVAL_NUM) %>%
    unique %>%
    bind_rows(partition_ts %>% dplyr::select(PATIENT_NUM) %>%
                dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)) %>%
                unique %>%
                dplyr::mutate(VARIABLE=paste0("ep_",j),
                              NVAL_NUM=1)) %>% 
    semi_join(data.frame(VARIABLE = model$feature_names,
                         stringsAsFactors = F),
              by="VARIABLE")
  
  #collect additional features required for training model
  x_add<-data.frame(VARIABLE = model$feature_names,
                    stringsAsFactors = F) %>%
    anti_join(data.frame(VARIABLE = as.character(unique(test$VARIABLE)),
                         stringsAsFactors = F),
              by="VARIABLE")
  
  #introduce psuedo-obs to align with training model
  if(nrow(x_add)>0){
    test %<>%
      bind_rows(data.frame(PATIENT_NUM = rep(0,nrow(x_add)),
                           VARIABLE = x_add$VARIABLE,
                           NVAL_NUM = 0,
                           stringsAsFactors = F))
  }

  test %<>%
    arrange(PATIENT_NUM) %>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM",
                          variable="VARIABLE",
                          val="NVAL_NUM")
  
  if(nrow(x_add)>0){
    test<-test[-1,]
  }
  
  #check alignment
  if(!all(colnames(test)==colnames(out$x_prep$X_ep))){
    warning("predictors don't align!")
  }
  
  #make prediction on validation set
  valid<-data.frame(PATIENT_NUM = row.names(test),
                    episode=j,
                    pred = predict(model,newdata=test),
                    stringsAsFactors = F)
  
  valid_stack2 %<>% bind_rows(valid)
  
  lapse<-Sys.time()-start
  cat("finish collecting predictions for year",j,"in",lapse,units(lapse),".\n")
}

partition_ts %<>%
  dplyr::mutate(PATIENT_NUM = as.character(PATIENT_NUM),
                real=DKD_IND_additive) %>%
  dplyr::select(PATIENT_NUM,episode,real,part_idx)

valid_stack2 %<>%
  left_join(partition_ts,by=c("PATIENT_NUM","episode")) %>%
  dplyr::filter(!((part_idx=="H") && (episode > 1))) %>%
  dplyr::mutate(real_ref=real) %>%
  group_by(PATIENT_NUM) %>%
  tidyr::fill(.,c("real_ref","part_idx"),.direction="up") %>%
  ungroup %>% dplyr::filter(!is.na(real_ref)) %>%
  replace_na(list(real=0)) %>%
  group_by(PATIENT_NUM,part_idx) %>%
  arrange(episode) %>%
  dplyr::mutate(pred2=cumsum(pred),
                pred3=cummax(pred),
                pred4=cumsum(1/pred),
                cum_row=1:n()) %>%
  ungroup %>%
  dplyr::mutate(pred4=cum_row/pred4,
                pred_wt2=pred2/cum_row) %>%
  group_by(episode) %>%
  dplyr::mutate(pred5=scales::rescale(pred2,c(0.1,0.99))) %>%
  dplyr::select(PATIENT_NUM,episode,pred,pred2,pred3,pred4,pred_wt2,
                real,part_idx) %>%
  dplyr::mutate(temporal=type_seq[i]) %>%
  semi_join(pat_episode %>% dplyr::select(PATIENT_NUM,episode) %>% 
              unique %>% dplyr::mutate(PATIENT_NUM=as.character(PATIENT_NUM)),
            by=c("PATIENT_NUM","episode"))

#quick AUC check
valid_stack2 %>%
  dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)) %>%
  group_by(episode,part_idx,temporal) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM*real))-1,
                   dkd_rate = round((length(unique(PATIENT_NUM*real))-1)/length(unique(PATIENT_NUM)),2),
                   low95=pROC::ci.auc(real,pred)[1],
                   auc=pROC::ci.auc(real,pred)[2],
                   up95=pROC::ci.auc(real,pred)[3]) %>%
  ungroup %>% arrange(part_idx, episode) %>%
  View

validation[[paste0("valid_stack",i)]]<-valid_stack2 %>% ungroup

####non-temporal####
i<-1
load(paste0("./data2/",time_iterv,"_",type_seq[i],"_gbm_model3.Rdata"))
load("./data/DKD_heron_facts_prep.Rdata")

fact_stack %<>% semi_join(partition_ts,by="PATIENT_NUM")

valid_stack1<-c()
for(j in 0:4){
  start<-Sys.time()
  
  model<-out$model_yr[[paste0(j,"yr_since_dm")]]
  feat_tr<-colnames(out$model_tr[[paste0(j,"yr_since_dm")]])
  
  #get predictors for test set
  test<-fact_stack %>% filter(yr_from_dm < j) %>%
    semi_join(partition_ts,by="PATIENT_NUM") %>%
    group_by(PATIENT_NUM, CONCEPT_CD) %>%
    top_n(n=-1,wt=day_to_end) %>% ungroup %>%
    dplyr::select(PATIENT_NUM, CONCEPT_CD, NVAL_NUM) %>%
    bind_rows(partition_ts %>% dplyr::select(PATIENT_NUM) %>%
                unique %>%
                dplyr::mutate(CONCEPT_CD=paste0("ep_",j),
                              NVAL_NUM=1)) %>%
    unique %>% semi_join(data.frame(CONCEPT_CD = feat_tr,
                                    stringsAsFactors = F),
                         by="CONCEPT_CD")
  
  #collect additional features required for training model
  x_add<-data.frame(CONCEPT_CD = feat_tr,
                    stringsAsFactors = F) %>%
    anti_join(data.frame(CONCEPT_CD = as.character(unique(test$CONCEPT_CD)),
                         stringsAsFactors = F),
              by="CONCEPT_CD")
  
  #introduce psuedo-obs to align with training model
  if(nrow(x_add) > 0){
    test %<>%
      bind_rows(data.frame(PATIENT_NUM = rep(0,nrow(x_add)),
                           CONCEPT_CD = x_add$CONCEPT_CD,
                           NVAL_NUM = 0,
                           stringsAsFactors = F))
  }

  test %<>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM",
                          variable="CONCEPT_CD",
                          val="NVAL_NUM")
  
  if(nrow(x_add) > 0){
    test<-test[-1,]
  }
  
  #check alignment
  if(!all(colnames(test)==feat_tr)){
    warning("predictors don't align!")
  }
  
  #make prediction on validation set
  valid<-data.frame(PATIENT_NUM = row.names(test),
                    episode=j,
                    pred = predict(model,newdata=test),
                    stringsAsFactors = F)
  
  valid_stack1 %<>% bind_rows(valid)
  
  lapse<-Sys.time()-start
  cat("finish collecting predictions for year",j,"in",lapse,units(lapse),".\n")
}

partition_ts %<>%
  dplyr::mutate(PATIENT_NUM = as.character(PATIENT_NUM),
                real=DKD_IND_additive) %>%
  dplyr::select(PATIENT_NUM,episode,real,part_idx)

valid_stack1 %<>%
  left_join(partition_ts,by=c("PATIENT_NUM","episode")) %>%
  dplyr::filter(!((part_idx=="H") && (episode > 1))) %>%
  dplyr::mutate(real_ref=real) %>%
  group_by(PATIENT_NUM) %>%
  tidyr::fill(.,c("real_ref","part_idx"),.direction="up") %>%
  ungroup %>% dplyr::filter(!is.na(real_ref)) %>%
  replace_na(list(real=0)) %>%
  group_by(PATIENT_NUM,part_idx) %>%
  arrange(episode) %>%
  dplyr::mutate(pred2=cumsum(pred),
                pred3=cummax(pred),
                pred4=cumsum(1/pred),
                cum_row=1:n()) %>%
  ungroup %>%
  dplyr::mutate(pred4=cum_row/pred4,
                pred_wt2=pred2/cum_row) %>%
  group_by(episode) %>%
  dplyr::mutate(pred5=scales::rescale(pred2,c(0.1,0.99))) %>%
  dplyr::select(PATIENT_NUM,episode,pred,pred2,pred3,pred4,pred_wt2,
                real,part_idx) %>%
  dplyr::mutate(temporal=type_seq[i]) %>%
  semi_join(pat_episode %>% dplyr::select(PATIENT_NUM,episode) %>% 
              unique %>% dplyr::mutate(PATIENT_NUM=as.character(PATIENT_NUM)),
            by=c("PATIENT_NUM","episode"))

#quick AUC check
valid_stack1 %>%
  dplyr::mutate(PATIENT_NUM=as.numeric(PATIENT_NUM)) %>%
  group_by(episode,part_idx,temporal) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM*real))-1,
                   dkd_rate = round((length(unique(PATIENT_NUM*real))-1)/length(unique(PATIENT_NUM)),2),
                   low95=pROC::ci.auc(real,pred)[1],
                   auc=pROC::ci.auc(real,pred)[2],
                   up95=pROC::ci.auc(real,pred)[3]) %>%
  ungroup %>% arrange(part_idx, episode) %>%
  View

validation[[paste0("valid_stack",i)]]<-valid_stack1 %>% ungroup

saveRDS(validation,file="./data2/validation.rda")

##================understand feature======================
load("./data2/1yr_non-temporal_gbm_model2.Rdata")
varimp1<-c()
for(i in 0:4){
  varimp1 %<>%
    bind_rows(xgb.importance(model=out$model_yr[[paste0(i,"yr_since_dm")]]) %>%
                dplyr::mutate(episode = paste0("yr_",i)) %>%
                dplyr::mutate(rank=rank(-Gain)) %>%
                dplyr::select(Feature,Gain,rank,episode))
}
varimp1 %<>%
  group_by(Feature) %>%
  dplyr::mutate(Gain_tot=sum(Gain),
                Gain_avg=mean(Gain),
                sel_ep=length(unique(episode)),
                rank_avg=mean(rank),
                rank_best=min(rank),
                ep_best=episode[which.min(rank)]) %>%
  ungroup %>%
  dplyr::select(Feature, Gain_tot, Gain_avg,
                sel_ep, rank_avg, rank_best, ep_best,
                Gain, episode) %>%
  spread(episode, Gain)

load("./data2/1yr_stack-temporal_gbm_model2.Rdata")
varimp2<-c()
for(i in 0:4){
  varimp2 %<>%
    bind_rows(xgb.importance(model=out$model_yr[[paste0(i,"yr_since_dm")]]) %>%
                dplyr::mutate(episode = paste0("yr_",i)) %>%
                dplyr::mutate(rank=rank(-Gain)) %>%
                dplyr::select(Feature,Gain,rank,episode))
}
varimp2 %<>%
  group_by(Feature) %>%
  dplyr::mutate(Gain_tot=sum(Gain),
                Gain_avg=mean(Gain),
                sel_ep=length(unique(episode)),
                rank_avg=mean(rank),
                rank_best=min(rank),
                ep_best=episode[which.min(rank)]) %>%
  ungroup %>%
  dplyr::select(Feature, Gain_tot, Gain_avg,
                sel_ep, rank_avg, rank_best, ep_best,
                Gain, episode) %>%
  spread(episode, Gain)


load("./data2/1yr_discrt-surv-temporal_gbm_model2.Rdata")
varimp3<-xgb.importance(model=out$model) %>%
  dplyr::mutate(rank=1:n())


load("./data2/1yr_dynamic-temporal_gbm_model2.Rdata")
varimp4<-c()
for(i in 0:4){
  varimp4 %<>%
    bind_rows(xgb.importance(model=out$model_yr[[paste0(i,"yr_since_dm_TRUE")]]) %>%
                dplyr::mutate(episode = paste0("yr_",i)) %>%
                dplyr::mutate(rank=rank(-Gain)) %>%
                dplyr::select(Feature,Gain,rank,episode))
}
varimp4 %<>%
  group_by(Feature) %>%
  dplyr::mutate(Gain_tot=sum(Gain),
                Gain_avg=mean(Gain),
                sel_ep=length(unique(episode)),
                rank_avg=mean(rank),
                rank_best=min(rank),
                ep_best=episode[which.min(rank)]) %>%
  ungroup %>%
  dplyr::select(Feature, Gain_tot, Gain_avg,
                sel_ep, rank_avg, rank_best, ep_best,
                Gain, episode) %>%
  spread(episode, Gain)

varimp<-list(varimp1=varimp1,
             varimp2=varimp2,
             varimp3=varimp3,
             varimp4=varimp4)

saveRDS(varimp,file="./data2/varimp.rda")



