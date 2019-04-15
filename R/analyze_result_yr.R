####Result Analysis ####

# setwd("~/proj_dkd/DKD_PM_temporal")

rm(list=ls()); gc()

source("./R/util.R")
require_libraries(c("tidyr",
                    "dplyr",
                    "magrittr",
                    "ggplot2",
                    "pROC",
                    "PRROC",
                    "ROCR",
                    "broom",
                    "purrr"))

##================prediction performance=====================
# temporal handling types
type_seq<-c(
  "non-temporal"
  ,"stack-temporal"
  ,"discrt-surv-temporal"
  ,"grp-temporal" #ignored for now
  ,"dynamic-temporal"
)

validation<-readRDS("./data2/validation.rda")
#separate between internal and external validation
valid_out<-c()
valid_inter<-c()

for(i in c(1,2,3,5)){
  dat<-validation[[paste0("valid_stack",i)]] %>%
    dplyr::mutate(PATIENT_NUM = as.character(PATIENT_NUM))
  
  # if(i %in% c(1,2)){
  #   dat<-dat %>%
  #     dplyr::select(PATIENT_NUM,episode,pred,real,part_idx,temporal)
  # }else{
  #   dat<-dat %>%
  #     dplyr::select(PATIENT_NUM,episode,pred_wt2,real,part_idx,temporal) %>%
  #     dplyr::rename(pred = pred_wt2)
  # }
  
  dat<-dat %>%
    dplyr::select(PATIENT_NUM,episode,pred_wt2,real,part_idx,temporal) %>%
    dplyr::rename(pred = pred_wt2)
  
  valid_out %<>% bind_rows(dat)
  
  if(i==1){
    valid_inter<<-dat 
  }else{
    dat %<>%
      semi_join(valid_inter,by="PATIENT_NUM")
  }
  
  valid_inter %<>% bind_rows(dat)
}

valid_inter %<>% unique %>%
  group_by(episode,temporal) %>%
  dplyr::mutate(pred=scales::rescale(pred,to=c(0,100))) %>%
  ungroup %>%
  dplyr::mutate(pred=round(pred)) %>%
  spread(temporal,pred) %>%
  na.omit()

valid_out %<>%
  group_by(PATIENT_NUM,episode,temporal) %>%
  top_n(n=1,wt=real) %>%
  ungroup


# collect performance measurement
# n_bin<-20
n_bin<-10

lp1<-unique(valid_out$temporal)
lp2<-unique(valid_out$episode)
# lp3<-unique(valid_out$part_idx)

perf_summ<-c()
perf_tbl<-c()
# calib_equal_cut<-c()
calib_equal_bin<-c()

for(i in lp1){
  for(j in lp2){
    dat_sub<-valid_out %>%
      dplyr::filter(temporal==i & episode==j)
    if(nrow(dat_sub)==0) next
    
    boots<-bootstrap(dat_sub,50)
    for(b in 1:50){
      dat_sub_b<-dat_sub %>% slice(attr(boots,"indices")[[b]])
      perf_summ %<>%
        bind_rows(get_perf_summ(pred=dat_sub_b$pred,
                                real=dat_sub_b$real)$perf_summ %>%
                    dplyr::mutate(boots=b,
                                  episode = j,
                                  temporal = i))
    }

    calib_equal_bin %<>%
      bind_rows(get_calibr(pred=dat_sub$pred,
                           real=dat_sub$real,
                           n_bin=n_bin) %>%
                  dplyr::mutate(episode = j,
                                temporal = i))
  }
}


# plot overall summary
perf_overall<-perf_summ %>% 
  filter(overall_meas %in% c("roauc",
                             "roauc_low",
                             "roauc_up",
                             "prauc1",
                             "mcc",
                             "acc",
                             "rec_sens",
                             "spec",
                             "ppv",
                             "npv")) %>%
  group_by(episode,temporal,overall_meas) %>%
  dplyr::summarize(val=median(meas_val),
                   low=quantile(meas_val,prob=0.025),
                   up=quantile(meas_val,prob=0.975)) %>%
  dplyr::mutate(overall_meas=recode(overall_meas,
                                    roauc="1. ROAUC",
                                    prauc1="2. PRAUC",
                                    acc="3. Accuracy",
                                    mcc="4. Matthews Correlation Coefficient",
                                    rec_sens="5. Sensitivity",
                                    spec="6. Specificity",
                                    ppv="7. Positive Predictive Value",
                                    npv="8. Negative Predictive Value"))


#significance of difference
pval<-valid_out %>% 
  dplyr::select(PATIENT_NUM,episode,temporal,pred) %>%
  group_by(episode,temporal) %>%
  dplyr::mutate(pred2=scales::rescale(pred,to=c(0,100))) %>%
  ungroup %>%
  dplyr::select(PATIENT_NUM,episode,temporal,pred2) %>%
  spread(temporal,pred2)

ttest<-c()
mtest<-c()
for(j in lp2){
  pval_j<-pval %>% filter(episode==j)
  ttest %<>%
    bind_rows(data.frame(`non-temporal`=t.test(pval_j$`dynamic-temporal`,
                                               pval_j$`non-temporal`,
                                               paired=T)[["p.value"]],
                         `stack-temporal`=t.test(pval_j$`dynamic-temporal`,
                                                 pval_j$`stack-temporal`,
                                                 paired=T)[["p.value"]],
                         `discrt-surv-temporal`=t.test(pval_j$`dynamic-temporal`,
                                                       pval_j$`discrt-surv-temporal`,
                                                       paired=T)[["p.value"]],
                         `dynamic-temporal`=t.test(pval_j$`dynamic-temporal`,
                                                   pval_j$`dynamic-temporal`,
                                                   paired=T)[["p.value"]],
                         episode=j))
  
  mtest %<>%
    bind_rows(data.frame(`non-temporal`=mcnemar.test(as.numeric(pval_j$`dynamic-temporal`>=75),
                                                     as.numeric(pval_j$`non-temporal`>=75))[["p.value"]],
                         `stack-temporal`=mcnemar.test(as.numeric(pval_j$`dynamic-temporal`>=75),
                                                       as.numeric(pval_j$`stack-temporal`>=75))[["p.value"]],
                         `discrt-surv-temporal`=mcnemar.test(as.numeric(pval_j$`dynamic-temporal`>=75),
                                                             as.numeric(pval_j$`discrt-surv-temporal`>=75))[["p.value"]],
                         `dynamic-temporal`=mcnemar.test(as.numeric(pval_j$`dynamic-temporal`>=75),
                                                         as.numeric(pval_j$`dynamic-temporal`>=75))[["p.value"]],
                         episode=j))
}

mtest %<>% 
  mutate(`dynamic.temporal`=pmax(`non.temporal`,
                                 `stack.temporal`,
                                 `discrt.surv.temporal`)) %>%
  dplyr::select(`dynamic.temporal`,episode) %>%
  gather(temporal,p_val,-episode) %>%
  mutate(temporal=gsub("\\.","-",temporal))


perf_overall %<>%
  left_join(mtest,by=c("episode","temporal")) %>%
  mutate(temporal=recode(temporal,
                         `non-temporal`="a.Most Recent Value",
                         `stack-temporal`="b.Stack Temporal",
                         `discrt-surv-temporal`="c.Discrete Survival",
                         `dynamic-temporal`="d.Landmark Boosting")) %>%
  replace_na(list(p_val=NA))


#Overall Performance comparison
ggplot(perf_overall,
       aes(x=episode,y=meas_val,
           linetype=temporal,shape=temporal))+
  geom_point(size=2)+geom_line()+
  # geom_point(data=perf_overall %>% filter(p_val<0.001),
  #            aes(x=episode,y=meas_val),
  #            color="red",shape=1,size=5)+
  # scale_shape(solid=F)+
  # scale_size(guide=FALSE)+
  # scale_color_discrete(guide=FALSE)+
  facet_wrap(~overall_meas,ncol=4,scales="free")+
  labs(x="Year since DM onset",
       y="Prediction Performance",
       linetype="Temporal Approach",
       shape="Temporal Approach")


#Calibration
calib_equal_bin %<>%
  mutate(temporal=recode(temporal,
                         `non-temporal`="a.Most Recent Value",
                         `stack-temporal`="b.Stack Temporal",
                         `discrt-surv-temporal`="c.Discrete Survival",
                         `dynamic-temporal`="d.Landmark Boosting")) %>%
  mutate(color=case_when(real_p > binCI_upper ~ round(real_p/binCI_upper,1)-1,
                         real_p < binCI_lower ~ round(real_p/binCI_lower,1)-1,
                         real_p >= binCI_lower | real_p <= binCI_upper ~ 0))

ggplot(calib_equal_bin,aes(x=episode,y=pred_bin,fill=color))+
  geom_tile()+
  scale_y_continuous(name="Cumulative proportion of observations (low risk -> high risk)",
                     breaks=seq_len(n_bin),
                     labels=seq(0.05,0.95,length.out = n_bin))+
  scale_fill_gradient2(breaks=c(-0.3,0,0.3),
                       labels=c("Underpredicted","Calibrated","Overpredicted"),
                       low = "darkgreen", mid = "white", high = "darkred")+
  # geom_rug(aes(size=expos),stat="identity",sides="rb")+
  scale_size_continuous(guide=F)+
  labs(x="# of Year since DM onset",fill="Calibration Bias")+
  facet_wrap(~temporal)

