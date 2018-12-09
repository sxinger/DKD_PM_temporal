#### Single Round of XGB ####
# for preliminary results build some confidence #
setwd("~/dkd/DKD_PM")

rm(list=ls()); gc()
source("../helper_functions.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"xgboost"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
                     ,"ggplot2"
))

## Load in cleaned-up data sets
resamples<-10
load("./data/DKD_heron_pats_prep.Rdata")
load("./data/DKD_heron_facts_prep.Rdata")

## Load partition index -- random
load(paste0("../DKD_FS/data/","random_sample",resamples,".Rdata"))
dat_sample<-dat_resample_rand[["resample1"]] %>%
  semi_join(pat_tbl %>% filter(year <= 2015 & year >= 2008),by="PATIENT_NUM")

#covert long df to wide sparse matrix (facts)
x_sparse_val<-fact_stack %>% 
  semi_join(dat_sample,by=c("PATIENT_NUM")) %>%
  long_to_sparse_matrix(.,id="PATIENT_NUM",
                        variable="CONCEPT_CD",
                        val="NVAL_NUM")

#covert long df to wide sparse matrix (patient-level facts)
x_sparse_pat<-pat_tbl %>% 
  semi_join(dat_sample,by=c("PATIENT_NUM")) %>%
  arrange(PATIENT_NUM)

y<-x_sparse_pat %>% dplyr::select(DKD_IND)
x_sparse_pat %<>%
  dplyr::select(-DKD_IND,-year) %>%
  gather(key,value,-PATIENT_NUM) %>%
  long_to_sparse_matrix(.,id="PATIENT_NUM",variable="key",val="value")


## Gradient Boosting (xgboost) + Logistic
#70-30
#separate training and testing sets and convert to xgb.DMatrix
trainx<-cbind(x_sparse_pat[(dat_sample$part73=="T"),which(colnames(x_sparse_pat)!="DKD_IND")],
              x_sparse_val[(dat_sample$part73=="T"),])
trainy<-as.vector(y[(dat_sample$part73=="T"),])

testx<-cbind(x_sparse_pat[(dat_sample$part73=="V"),which(colnames(x_sparse_pat)!="DKD_IND")],
             x_sparse_val[(dat_sample$part73=="V"),])
testy<-as.vector(y[(dat_sample$part73=="V"),])

dtrain<-xgb.DMatrix(data=trainx,label=trainy)
dtest<-xgb.DMatrix(data=testx,label=testy)


#tune hyper-parameters
eta_rg<-c(0.1,0.05,0.01)
maxdep_rg<-c(2,4,6,8,10)
gamma_rg<-c(1,2)
colsample_rg<-c(0.8)
mcw_rg<-c(1)
eval_metric<-"auc"
# eval_metric<-"poisson-nloglik"
nrounds<-2000
objective<-"binary:logistic"
# objective<-"count:poisson"
verb<-TRUE

grid_params<-expand.grid(max.depth=maxdep_rg,
                         eta=eta_rg,
                         min_child_weight=mcw_rg,
                         subsample=0.8,
                         colsample_bytree=colsample_rg,
                         gamma=gamma_rg,
                         verb=T)
bst_grid<-c()
metric_name<-paste0("test_", eval_metric,"_mean")
metric_sd_name<-paste0("test_", eval_metric,"_std")

for(i in seq_len(dim(grid_params)[1])){
  start_i<-Sys.time()
  param<-as.list(grid_params[i,])
  param$eval_metric<-eval_metric
  
  bst <- xgb.cv(param,
                dtrain,
                nrounds=nrounds,
                nfold = 5,
                objective = objective,
                metrics = eval_metric,
                maximize = TRUE,
                early_stopping_rounds = 200,
                print_every_n = 100)
  bst_grid<-rbind(bst_grid, cbind(grid_params[i,],
                                  metric=max(bst$evaluation_log[[metric_name]]),
                                  steps=which(bst$evaluation_log[[metric_name]]==max(bst$evaluation_log[[metric_name]]))[1]))
  
  if(verb){
    cat('finished train case:',paste(c(colnames(grid_params),eval_metric),param,collapse="_"),
        'in',Sys.time()-start_i,units(Sys.time()-start_i))
    start_i<-Sys.time()
  }
}
hyper_param<-bst_grid[which.max(bst_grid$metric),]

#re-train using the tuned hyper-parameters
watchlist<-list(train=dtrain, test=dtest)
xgb_tune<-xgb.train(data=dtrain,
                    max_depth=hyper_param$max.depth,
                    maximize = TRUE,
                    # max_depth=10,
                    eta=hyper_param$eta,
                    # eta=0.01,
                    nrounds=hyper_param$steps,
                    # nrounds=2000,
                    # early_stopping_rounds = 200,
                    watchlist=watchlist,
                    # eval_metric="error",eval_metric="logloss",
                    eval_metric="auc",
                    objective="binary:logistic")

#overall prediction
pred_real<-data.frame(pred=as.numeric(predict(xgb_tune,dtest)),
                      real=testy)
pROC::auc(roc(pred_real$real,pred_real$pred))
pROC::ci.auc(roc(pred_real$real,pred_real$pred))
pROC::coords(roc(pred_real$real,pred_real$pred),x="best")

#decile table
decile_table<-pred_real %>%
  dplyr::arrange(pred) %>%
  mutate(risk_grp=cut(pred,breaks=quantile(pred,probs=0:10/10),
                      include.lowest=T,labels=F)) %>%
  group_by(risk_grp) %>%
  dplyr::summarize(lower=min(pred),med=median(pred),upper=max(pred),expos=n(),
                   DKD_pred = sum(pred),DKD=sum(real), DKD_rt=mean(real,na.rm=T),
                   nonDKD_pred = n()-sum(pred), non_DKD=n()-sum(real),non_DKD_rt=1-mean(real,na.rm=T)) %>%
  mutate(lift = round(DKD_rt/mean(pred_real$real),2))

#Hosmer Lemeshow
lrisk_sep<-c()
hrisk_sep<-c()
hl_score<-c()

for(grp in seq(10,100,by=10)){
  temp<-pred_real %>%
    dplyr::arrange(pred) %>%
    mutate(risk_grp=cut(pred,breaks=quantile(pred,probs=0:grp/grp),
                        include.lowest=T,labels=F)) %>%
    group_by(risk_grp) %>%
    dplyr::summarize(DKD_rt=mean(real),
                     Chi_sq1=(sum(real)-sum(pred))^2/(sum(pred)),
                     Chi_sq0=(sum(pred)-sum(real))^2/(n()-sum(pred)))
  DKD_rt<-unlist(temp$DKD_rt)
  lrisk_sep<-c(lrisk_sep,sum(DKD_rt[1]/0.5))
  hrisk_sep<-c(hrisk_sep,sum(DKD_rt[grp]/0.5))
  
  temp %<>%
    gather(key, value,-risk_grp,-DKD_rt) %>%
    dplyr::summarize(Chi_sq=sum(value),df=n()) %>%
    mutate(p_val=pchisq(Chi_sq,df))
  hl_score<-c(hl_score,temp$p_val)
}
out<-data.frame(risk_grp=seq(10,100,by=10),
                low_risk_multiple=lrisk_sep,
                high_risk_multiple=hrisk_sep,
                hl_score=hl_score)


#get variable importance matrics
#temp variable dictionary
feat_dict<-fact_stack %>% 
  group_by(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, C_NAME) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM))) %>%
  bind_rows(pat_tbl %>%
              dplyr::select(-DKD_IND) %>%
              gather(CONCEPT_CD,value,-PATIENT_NUM) %>%
              dplyr::filter(value!=0) %>%
              group_by(CONCEPT_CD) %>%
              dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM))) %>%
              mutate(VARIABLE_CATEG="DEMOGRAPHICS",
                     C_VISUAL_PATH="patient_dimension",
                     C_NAME=CONCEPT_CD) %>%
              dplyr::select(VARIABLE_CATEG,C_VISUAL_PATH,CONCEPT_CD, C_NAME,pat_cnt))

varimp_xgb_tune_val<-xgb.importance(colnames(trainx),model=xgb_tune) %>%
  dplyr::filter(Gain > 0) %>%
  left_join(feat_dict, by=c("Feature"="CONCEPT_CD")) %>% unique #decode


imp_var_dist<-varimp_xgb_tune_val %>% 
  mutate(rank=1:nrow(varimp_xgb_tune_val),
         cum_gain = cumsum(Gain)) %>%
  mutate(gain_80perc = ifelse(cum_gain<=0.8, Gain, 0)) %>%
  group_by(VARIABLE_CATEG) %>%
  dplyr::summarize(concept_cnt=n(),
                   best_rank=min(rank),
                   mid_rank=median(rank),
                   cumulant_gain=sum(Gain),
                   cumulant_gain_80perc=sum(gain_80perc),
                   top100_cnt=sum((rank<=100)))


#write up and save the results
output_date<-format(Sys.Date(),"%m%d%Y")
write.csv(varimp_xgb_tune_val,file=paste0("./data/single_xgb_fs_",output_date,".csv"),row.names=F)


#Profiling
pred_real %<>%
  dplyr::mutate(PATIENT_NUM = pat_tbl$PATIENT_NUM[(dat_sample$part73=="V")])

pat_tbl_long<-pat_tbl %>%
  gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM) %>%
  mutate(C_NAME=CONCEPT_CD) %>% 
  select(PATIENT_NUM,CONCEPT_CD,C_NAME,NVAL_NUM)

top<-pred_real %>%
  filter(pred>=0.8) %>%
  left_join(fact_stack %>% select(PATIENT_NUM,CONCEPT_CD,C_NAME,NVAL_NUM) %>%
              bind_rows(pat_tbl_long),
            by="PATIENT_NUM") %>%
  semi_join(varimp_xgb_tune_val %>% dplyr::slice(1:100),
            by=c("CONCEPT_CD"="Feature")) %>%
  mutate(risk_grp="top")
  
bottom<-pred_real %>%
  filter(pred<=0.06) %>%
  left_join(fact_stack %>% select(PATIENT_NUM,CONCEPT_CD,C_NAME,NVAL_NUM) %>%
              bind_rows(pat_tbl_long),
            by="PATIENT_NUM") %>%
  semi_join(varimp_xgb_tune_val %>% dplyr::slice(1:20),
            by=c("CONCEPT_CD"="Feature")) %>%
  mutate(risk_grp="bottom")



top_bottom<-top %>% bind_rows(bottom) %>%
  group_by(risk_grp) %>%
  dplyr::mutate(risk_grp_size=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  group_by(C_NAME,risk_grp) %>%
  dplyr::summarize(risk_grp_size=risk_grp_size[1],
                   expos=length(unique(PATIENT_NUM)),
                   mean=mean(NVAL_NUM,na.rm=T),
                   sd=sd(NVAL_NUM,na.rm=T),
                   sum=sum(NVAL_NUM)) %>%
  dplyr::mutate(mean=ifelse(mean==1 & sd==0,round(sum/risk_grp_size,3),mean),
                sd=ifelse(sd==0,sqrt((sum/risk_grp_size)*(1-(sum/risk_grp_size))),sd)) %>%
  select(C_NAME,risk_grp,expos,mean,sd)




## Partial Plots for the Top ith variables
#partial effects
require_libraries("pdp")
feature_rk<-2
pred_var<-varimp_xgb_tune_val$Feature[feature_rk]
var_chk<-trainx[which(trainx[,pred_var]>0, arr.ind=T),pred_var]

pred_qtl<-quantile(var_chk[order(var_chk)],probs=0:10/10) #Age, BMI, SBP, weight, BUN, DBP
# pred_qtl<-quantile(var_chk[order(var_chk)],probs=0:6/6) #bilirubin, MPV
# pred_qtl<-quantile(var_chk[order(var_chk)],probs=c(0,0.5,0.8,0.9,1)) #infection counts
pred_qtl<-c(floor(pred_qtl[1]),round(pred_qtl[-c(1,length(pred_qtl))],2),ceiling(pred_qtl[length(pred_qtl)]))

#numerical
expos<-data.frame(pred_var = var_chk) %>%
  dplyr::arrange(pred_var) %>%
  mutate(var_grp=cut(pred_var,breaks=pred_qtl,right=F,include.lowest=T)) %>% 
  group_by(var_grp) %>% 
  dplyr::summarize(group=row_number(),
            midpt=median(pred_var),
            exposure=n())

#categorical
expos<-data.frame(pred_var = trainx[,pred_var]) %>%
  mutate(var_grp=as.factor(pred_var)) %>% 
  group_by(var_grp) %>% 
  dplyr::summarize(midpt=median(pred_var),
            exposure=n())

get_peff2<-partial_effect( fit=xgb_tune
                           ,pred_var=pred_var
                           ,pred_grid=expos$midpt
                           ,x=trainx
                           ,nrs=30)

get_peff2 <- get_peff2 %>%
  mutate(pred_var = get(pred_var)) %>%
  arrange(pred_var,yhat) %>%
  group_by(pred_var) %>%
  dplyr::summarize(lower95=quantile(yhat,prob=0.05,na.rm=T),
            mean=mean(yhat,na.rm=T),
            upper95=quantile(yhat,prob=0.95,na.rm=T)) %>%
  left_join(expos,by=c("pred_var"="midpt"))


pplot<-ggplot2::ggplot(get_peff2,aes(x=pred_var,y=mean))+
  ggplot2::geom_line(aes(y=mean))+
  ggplot2::geom_point(aes(y=mean,size=exposure),colour="red")+
  ggplot2::geom_ribbon(aes(ymin=lower95,
                           ymax=upper95),alpha=0.2)+
  labs(x=varimp_xgb_tune_val$C_NAME[feature_rk],
       y="Predictted Probability")

if(nrow(get_peff2)>2){
  pplot<-pplot+
    scale_x_discrete(name=varimp_xgb_tune_val$C_NAME[feature_rk], 
                     limits=c(pred_qtl,expos$midpt)[order(c(pred_qtl,expos$midpt))])+
    theme(axis.text.x = element_text(face="bold", color="#993333",size=12, angle=45))
}

pplot
