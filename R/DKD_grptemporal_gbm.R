setwd("~/dkd/DKD_PM")

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

# Load data
load("./data/X_long.Rdata")
load("./data/pat_episode.Rdata")

#==== evaluations
# time_iterv<-"3mth"
# time_iterv<-"6mth"
time_iterv<-"1yr"

# ep_unit<-90
# ep_unit<-182.5
ep_unit<-365.25
pat_episode<-pat_tbl %>%
  dplyr::mutate(episode = floor(as.numeric(DAY_SINCE)/ep_unit))

partition<-pat_episode %>%
  semi_join(X_long, by = "PATIENT_NUM") %>%
  dplyr::select(PATIENT_NUM,part_idx) %>%
  unique

# nfold for cv
nfold<-5 

# 
type<-"grp-temporal"

#
start_k<-Sys.time()
cat("...temporal data integrated by:",type,'\n')

bm<-c()
bm_nm<-c()
#==build training set
start_k_i<-Sys.time()

Xy_all<-get_grptemporal(pat_episode,X_long)

dtrain<-xgb.DMatrix(data=Xy_all$X_ep,
                    label=Xy_all$y_ep$DKD_IND_additive[!(partition$part_idx %in% c("H","5"))])

# dtest<-xgb.DMatrix(data=Xy_all$X_ep[(partition$part_idx %in% c("H","5")),],
#                    label=Xy_all$y_ep$DKD_IND_additive[(partition$part_idx %in% c("H","5"))])

lapse_k_i<-Sys.time()-start_k_i
cat("......finish handling temporal data in",lapse_k_i,units(lapse_k_i),"\n")
bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
bm_nm<-c(bm_nm,"build temporal traning set")

# #==get indices for k folds
# folds<-list()
# train_dt<-train_pt %>%
#   dplyr::select(PATIENT_NUM,part_idx) %>% unique %>%
#   dplyr::mutate(row_idx = 1:n())
# 
# for(fd in seq_len(nfold)){
#   fd_df<-train_dt %>% filter(part_idx==as.character(fd)) %>%
#     dplyr::select(row_idx)
#   folds[[fd]]<-fd_df$row_idx
# }
#==tune hyperparameter
#hyper-parameter grid for xgboost
eval_metric<-"auc"
objective<-"binary:logistic"
grid_params_tree<-expand.grid(
  # max_depth=c(4,6,8,10),
  max_depth=8,
  # eta=c(0.05,0.02,0.01),
  eta=0.01,
  min_child_weight=1,
  subsample=0.8,
  colsample_bytree=0.8, 
  gamma=1
)

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
  # param$scale_pos_weight=mean(train$y_train$DKD_IND_additive) #inbalance sampling
  param$scale_pos_weight=1 #balance sampling
  
  bst <- xgb.cv(param,
                dtrain,
                objective = objective,
                metrics = eval_metric,
                maximize = TRUE,
                nrounds=2000,
                nfold = 5,
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
valid_cv<-data.frame(PATIENT_NUM = row.names(Xy_all$X_ep),
                     valid_type = 'cv',
                     pred = bst_grid_cv[,which.max(bst_grid$metric)],
                     real = getinfo(dtrain,"label"),
                     method = type,
                     stringsAsFactors = F)

lapse_k_i<-Sys.time()-start_k_i
cat("......finish tunning hyperparameter in",lapse_k_i,units(lapse_k_i),"\n")
bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
bm_nm<-c(bm_nm,"tune hyperparameters")

#==re-train model for external validation
start_k_i<-Sys.time()

xgb_tune<-xgb.train(data=dtrain,
                    max_depth=hyper_param$max_depth,
                    maximize = TRUE,
                    eta=hyper_param$eta,
                    nrounds=hyper_param$steps,
                    eval_metric="auc",
                    objective="binary:logistic",
                    print_every_n = 100)

lapse_k_i<-Sys.time()-start_k_i
cat("......finish re-training in",lapse_k_i,units(lapse_k_i),"\n")
bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
bm_nm<-c(bm_nm,"re-training")

lapse_k<-Sys.time()-start_k
cat("...finish evaluating temporal integraion method:",type,"in",lapse_k,units(lapse_k),"\n")
bm<-c(bm,paste0(round(lapse_k_i,4),units(lapse_k_i)))
bm_nm<-c(bm_nm,"end method loop")


#==save results
bm<-data.frame(step=bm_nm,time=bm)
out<-list(x_prep = Xy_all,
          pred_cv = valid_cv,
          model = xgb_tune,
          bm = bm)

save(out,file=paste0("./data/",time_iterv,"_",type,"_gbm_model.Rdata"))

rm(list=ls())
gc()
