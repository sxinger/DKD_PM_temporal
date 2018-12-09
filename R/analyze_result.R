####Result Analysis ####
setwd("C:/Users/xsong/Desktop/Projects/DKD/Sci_Report")

rm(list=ls()); gc()

# install.packages("PRROC")
# install.packages("ROCR")
# install.packages("rlang")
# install.packages("xgboost")
# install.packages("ResourceSelection")

library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(pROC)
library(PRROC)
library(ROCR)
library(xgboost)
library(ResourceSelection)

# temporal handling types
type_seq<-c(
  "non-temporal"
 ,"stack-temporal"
 ,"discrt-surv-temporal"
 ,"grp-temporal"
 ,"dynamic-temporal"
)


#separate between internal and external validation
valid_out<-c()
valid_inter<-c()

for(i in c(1,2,3,5)){
  load(paste0("./data/valid_stack",i,".Rdata"))
  dat<-get(paste0("valid_stack",i)) %>%
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
    dplyr::select(PATIENT_NUM,episode,pred,real,part_idx,temporal)
  
  # dat<-dat %>%
  #   dplyr::select(PATIENT_NUM,episode,pred_wt2,real,part_idx,temporal) %>%
  #   dplyr::rename(pred = pred_wt2)

  valid_out %<>% bind_rows(dat)
  
  if(i==1){
    valid_inter<<-dat 
  }else{
    dat %<>%
      semi_join(valid_inter,by="PATIENT_NUM")
  }
  
  valid_inter %<>% bind_rows(dat)
}

valid_inter %<>% unique %>% unique %>%
  group_by(episode,temporal) %>%
  dplyr::mutate(pred=scales::rescale(pred,to=c(0,100))) %>%
  ungroup %>% spread(temporal,pred) %>%
  na.omit()


# collect performance measurement
n_bin<-20
# n_bin<-100

lp1<-unique(valid_out$temporal)
# lp2<-unique(valid_out$episode)
# lp3<-unique(valid_out$part_idx)

perf_summ<-c()
perf_tbl<-c()
calib_equal_cut<-c()
calib_equal_bin<-c()

for(i in lp1){
  dat_sub<-valid_out %>% ungroup %>%
    dplyr::filter(temporal==i) %>%
    group_by(PATIENT_NUM) %>%
    top_n(n=1,wt=episode)
  
  if(nrow(dat_sub)==0) next
  
  # various performace table
  pred<-ROCR::prediction(dat_sub$pred,
                         dat_sub$real)
  
  prc<-performance(pred,"prec","rec")
  roc<-performance(pred,"sens","spec")
  nppv<-performance(pred,"ppv","npv")
  
  pcfall<-performance(pred,"pcfall")
  acc<-performance(pred,"acc")
  fscore<-performance(pred,"f")
  # mxe<-performance(pred,"mxe")
  # rmse<-performance(pred,rmse)
  mcc<-performance(pred,"phi")

  perf_at<-data.frame(cutoff=prc@alpha.values[[1]],
                      prec=prc@y.values[[1]],
                      rec_sens=prc@x.values[[1]],
                      stringsAsFactors = F) %>% 
    arrange(cutoff) %>%
    left_join(data.frame(cutoff=nppv@alpha.values[[1]],
                         ppv=nppv@y.values[[1]],
                         npv=nppv@x.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(prec_rec_dist=abs(prec-rec_sens)) %>%
    left_join(data.frame(cutoff=fscore@x.values[[1]],
                         fscore=fscore@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=roc@alpha.values[[1]],
                         spec=roc@x.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(Euclid_meas=sqrt((1-rec_sens)^2+(0-(1-spec))^2),
                  Youden_meas=rec_sens+spec-1) %>%
    left_join(data.frame(cutoff=pcfall@x.values[[1]],
                         pcfall=pcfall@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=acc@x.values[[1]],
                         acc=acc@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    left_join(data.frame(cutoff=mcc@x.values[[1]],
                         mcc=mcc@y.values[[1]],
                         stringsAsFactors = F),
              by="cutoff") %>%
    dplyr::mutate(temporal = i) %>%
    filter(prec > 0 & rec_sens > 0 & spec > 0)
  
  perf_tbl %<>% bind_rows(perf_at)
  
  # performance summary
  lab1<-dat_sub$pred[dat_sub$real==1]
  lab0<-dat_sub$pred[dat_sub$real==0]
  pr<-pr.curve(scores.class0 = lab1,
               scores.class1 = lab0,curve=F)
  roc_ci<-pROC::ci.auc(dat_sub$real,dat_sub$pred)
  hl<-hoslem.test(dat_sub$pred,dat_sub$real)
  
  perf_summ<-perf_summ %>%
    bind_rows(data.frame(overall_meas=c("roauc_low",
                                        "roauc",
                                        "roauc_up",
                                        "opt_thresh",
                                        "opt_sens",
                                        "opt_spec",
                                        "prauc1",
                                        "prauc2",
                                        "opt_prec",
                                        "opt_rec",
                                        "opt_fscore",
                                        "calib_chi",
                                        "calib_p",
                                        "size"),
                         meas_val=c(roc_ci[[1]],
                                    roc_ci[[2]],
                                    roc_ci[[3]],
                                    perf_at$cutoff[which.min(perf_at$Euclid_meas)],
                                    perf_at$rec_sens[which.min(perf_at$Euclid_meas)],
                                    perf_at$spec[which.min(perf_at$Euclid_meas)],
                                    pr$auc.integral,
                                    pr$auc.davis.goadrich,
                                    perf_at$prec[which.min(perf_at$prec_rec_dist)],
                                    perf_at$rec_sens[which.min(perf_at$prec_rec_dist)],
                                    perf_at$fscore[which.min(perf_at$prec_rec_dist)],
                                    hl[["statistic"]],
                                    hl[["p.value"]],
                                    length(unique(dat_sub$PATIENT_NUM))),
                         stringsAsFactors = F) %>%
                dplyr::mutate(temporal = i))
  
  # calibration 1
  calib1<-data.frame(pred=dat_sub$pred,
                     real=dat_sub$real) %>%
    arrange(pred) %>%
    dplyr::mutate(pred=scales::rescale(pred,to=c(0,1))) %>%
    dplyr::mutate(pred_bin = cut(pred,
                                 breaks=seq(0,1,length.out = n_bin+1),
                                 include.lowest=T,
                                 labels=F)) %>%
    ungroup %>% group_by(pred_bin) %>%
    dplyr::summarize(expos=n(),
                     bin_lower=min(pred),
                     bin_upper=max(pred),
                     bin_mid=median(pred),
                     real_agg = sum(real),
                     pred_p = mean(pred),
                     pred_p_sd = sd(pred)) %>%
    dplyr::mutate(real_p=real_agg/expos) %>%
    dplyr::mutate(binCI_lower = pmax(0,pred_p-1.96*sqrt(real_p*(1-real_p)/expos)),
                  binCI_upper = pred_p+1.96*sqrt(real_p*(1-real_p)/expos)) %>%
    dplyr::mutate(temporal = i)
  
  calib_equal_cut %<>% bind_rows(calib1)
  
  # calibration 2
  calib2<-data.frame(pred=dat_sub$pred,
                     real=dat_sub$real) %>%
    arrange(pred) %>%
    dplyr::mutate(pred_bin = cut(pred,
                                 breaks=unique(quantile(pred,0:(n_bin)/(n_bin))),
                                 include.lowest=T,
                                 labels=F)) %>%
    ungroup %>% group_by(pred_bin) %>%
    dplyr::summarize(expos=n(),
                     bin_lower=min(pred),
                     bin_upper=max(pred),
                     bin_mid=median(pred),
                     real_agg = sum(real),
                     pred_p = mean(pred),
                     pred_p_sd = sd(pred)) %>%
    dplyr::mutate(real_p=real_agg/expos) %>%
    dplyr::mutate(binCI_lower = pmax(0,pred_p-1.96*sqrt(real_p*(1-real_p)/expos)),
                  binCI_upper = pred_p+1.96*sqrt(real_p*(1-real_p)/expos)) %>%
    dplyr::mutate(temporal = i)
  
  calib_equal_bin %<>% bind_rows(calib2)
}

perf_summ %<>%
  bind_rows(perf_tbl %>% group_by(temporal) %>%
              dplyr::summarize(prec=mean(prec),
                               sens=mean(rec_sens),
                               spec=mean(spec),
                               ppv=mean(ppv),
                               npv=mean(npv),
                               acc=mean(acc),
                               fscore=mean(fscore),
                               mcc=mean(mcc)) %>%
              gather(overall_meas,meas_val,-temporal) %>%
              dplyr::select(overall_meas,meas_val,temporal) %>%
              ungroup)


# plot overall summary
perf_overall<-perf_summ %>% 
  filter(overall_meas %in% c("roauc",
                             "prauc1",
                             "mcc",
                             "acc",
                             "sens",
                             "spec",
                             "ppv",
                             "npv",
                             "calib_chi",
                             "calib_p")) %>%
  dplyr::mutate(overall_meas=recode(overall_meas,
                                    roauc="01. ROAUC",
                                    prauc1="02. PRAUC",
                                    acc="03. Accuracy",
                                    calib_chi="04. Calibration",
                                    calib_p="05. Calibration(p-value)",
                                    mcc="06. Matthews Correlation Coefficient",
                                    sens="07. Average Sensitivity",
                                    spec="08. Average Specificity",
                                    ppv="09. Average Positive Predictive Value",
                                    npv="10. Average Negative Predictive Value"),
                temporal=recode(temporal,
                                `non-temporal`="a.Most Recent Value",
                                `stack-temporal`="b.Stack Temporal",
                                `discrt-surv-temporal`="c.Discrete Survival",
                                `dynamic-temporal`="d.Self Boosting")) %>%
  spread(temporal,meas_val)


# plot complete ROC and PRC and mark optimal points
perf_tbl %<>%
  mutate(temporal=recode(temporal,
                         `non-temporal`="a.Most Recent Value",
                         `stack-temporal`="b.Stack Temporal",
                         `discrt-surv-temporal`="c.Discrete Survival",
                         `dynamic-temporal`="d.Self Boosting"))

  

perf_roc<-perf_tbl %>% 
  dplyr::select(cutoff,rec_sens,spec,temporal) %>%
  dplyr::rename(y=rec_sens) %>%
  mutate(x=1-spec) %>%
  dplyr::select(cutoff,x,y,temporal)

perf_roc_opt<-perf_tbl %>%
  group_by(temporal) %>%
  filter(Euclid_meas==min(Euclid_meas)) %>%
  mutate(x=1-spec,y=rec_sens,label=round(cutoff,2)) %>%
  dplyr::select(temporal,x,y,label) %>%
  ungroup %>%
  dplyr::mutate(vjust=c(1,0,-1,-2))
  
p1<-ggplot(perf_roc,aes(x=x,y=y,linetype=temporal,color=temporal))+
  geom_line(size=1)+
  geom_point(data=perf_roc_opt,aes(x=x,y=y))+
  geom_text(data=perf_roc_opt,aes(x=x,y=y,label=label,vjust=vjust))+
  geom_abline(slope=1,intercept=0)+
  labs(x="1-Specificity",
       y="Sensitivity",
       color="Temporal Handler",
       linetype="Temporal Handler",
       title="Reciever Operation Curve")

perf_prc<-perf_tbl %>% 
  dplyr::select(cutoff,rec_sens,prec,temporal) %>%
  dplyr::rename(x=rec_sens) %>%
  mutate(y=prec) %>%
  dplyr::select(cutoff,x,y,temporal)

perf_prc_opt<-perf_tbl %>%
  group_by(temporal) %>%
  filter(prec_rec_dist==min(prec_rec_dist)) %>%
  mutate(x=rec_sens,y=prec,label=round(cutoff,2)) %>%
  dplyr::select(temporal,x,y,label) %>%
  ungroup %>%
  dplyr::mutate(vjust=c(1,0,-1,-2))

p2<-ggplot(perf_prc,aes(x=x,y=y,linetype=temporal,color=temporal))+
  geom_line(size=1)+
  geom_point(data=perf_prc_opt,aes(x=x,y=y))+
  geom_text(data=perf_prc_opt,aes(x=x,y=y,label=label,vjust=vjust))+
  ylim(0,1)+
  labs(x="Recall",
       y="Precision",
       color="Temporal Handler",
       linetype="Temporal Handler",
       title="Precision Recall Curve")


# install.packages("ggpubr")
library(ggpubr)
ggarrange(p1,p2,ncol=2,common.legend = TRUE, legend="right")


# plot calibrations
calib_equal_bin %<>%
  mutate(temporal=recode(temporal,
                         `non-temporal`="a.Most Recent Value",
                         `stack-temporal`="b.Stack Temporal",
                         `discrt-surv-temporal`="c.Discrete Survival",
                         `dynamic-temporal`="d.Self Boosting")) %>%
  mutate(color=case_when(real_p > binCI_upper ~ round(real_p/binCI_upper,1),
                         real_p < binCI_lower ~ round(real_p/binCI_lower,1),
                         real_p >= binCI_lower | real_p <= binCI_upper ~ 1))

p3<-ggplot(calib_equal_bin,aes(x=real_p,y=pred_p,color=temporal))+
  geom_point()+ geom_smooth(method="loess",alpha=0.2)+
  geom_abline(slope=1,intercept=0)+
  xlim(0,1)+ylim(0,1)+
  labs(x="Actual Probability",
       y="Predictted Probability",
       title="Calibration Plot",
       color="Temporal Handler")
p3

p4<-ggplot(calib_equal_bin,aes(x=as.factor(pred_bin),y=real_p,fill=temporal))+
  geom_bar(stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=binCI_lower,ymax=binCI_upper),
  #               width=0.2,position=position_dodge(0.9))+
  labs(x="Predictted Risk Group",
       y="Acutal Probability",
       title="Decile Table",
       fill="Temporal Handler")
p4


# plot calibration
calib_equal_cut %<>%
  mutate(temporal=recode(temporal,
                         `non-temporal`="a.Most Recent Value",
                         `stack-temporal`="b.Stack Temporal",
                         `discrt-surv-temporal`="c.Discrete Survival",
                         `dynamic-temporal`="d.Self Boosting")) %>%
  mutate(color=case_when(real_p > binCI_upper ~ round(real_p/binCI_upper,1)-1,
                         real_p < binCI_lower ~ round(real_p/binCI_lower,1)-1,
                         real_p >= binCI_lower | real_p <= binCI_upper ~ 0))

ggplot(calib_equal_cut,aes(x=temporal,y=pred_bin,fill=color))+
  geom_tile()+
  scale_y_continuous(name="Predictted Probability",
                     breaks=seq_len(n_bin),
                     labels=seq(0.05,0.95,length.out = n_bin))+
  scale_fill_gradient2(breaks=c(-0.5,0,0.5),
                       labels=c("Underpredicted","Calibrated","Overpredicted"),
                       low = "darkgreen", mid = "white", high = "darkred")+
  geom_rug(aes(size=expos),stat="identity",sides="rb")
