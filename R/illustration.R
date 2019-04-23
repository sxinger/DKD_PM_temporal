## illustrations ##
# setwd("~/DKD_PM_temporal")

rm(list=ls())

source("./R/util.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
                     ,"ggplot2"
                     ,"ggthemes"
))

fact_stack<-readRDS("./data2/DKD_heron_facts_prep.rda") %>%
  anti_join(readRDS("./data2/pat_T1DM.rda"),by="PATIENT_NUM")

pat_tbl<-readRDS("./data2/pat_episode2.rda") %>%
  anti_join(readRDS("./data2/pat_T1DM.rda"),by="PATIENT_NUM") %>%
  dplyr::mutate(tr_ts=ifelse(part_idx==5,"ts","tr"))

pat_tbl %>%
  group_by(tr_ts) %>%
  dplyr::summarise(pat_cnt=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  View

ep_unit<-365.25
pat_tbl %>%
  dplyr::mutate(episode = floor(as.numeric(DAY_SINCE)/ep_unit)) %>%
  group_by(tr_ts,episode) %>%
  dplyr::summarise(pat_cnt=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  View


#dkd rates
N<-length(unique(pat_tbl$PATIENT_NUM))
dkd_rt<-pat_tbl %>%
  dplyr::select(PATIENT_NUM,DKD_IND_additive) %>%
  unique %>%
  group_by(PATIENT_NUM) %>%
  top_n(n=1L,wt=DKD_IND_additive) %>%
  ungroup %>%
  group_by(DKD_IND_additive) %>%
  dplyr::summarise(pat_cnt=n(),
                   pat_prop=round(pat_cnt/N,4))

cd_dist<-fact_stack %>%
  group_by(VARIABLE_CATEG) %>%
  dplyr::summarize(fact_cnt=n(),
                   pat_cnt=length(unique(PATIENT_NUM)),
                   cd_cnt=length(unique(CONCEPT_CD)),
                   pat_prop=round(length(unique(PATIENT_NUM))/N,2)) %>%
  ungroup

#overall demographics
demo<-fact_stack %>% filter(VARIABLE_CATEG=="demographics") %>%
  # filter(NVAL_NUM==1|CONCEPT_CD=="AGE_at_DM") %>%
  left_join(pat_tbl %>% group_by(PATIENT_NUM) %>%
              top_n(n=1L,wt=DAY_SINCE) %>%
              top_n(n=1L,wt=DKD_IND_additive) %>%
              dplyr::select(PATIENT_NUM,DKD_IND_additive) %>%
              unique,
            by="PATIENT_NUM")

#age in numeric
demo %>%
  filter(CONCEPT_CD=="AGE_at_DM") %>%
  dplyr::select(NVAL_NUM,DKD_IND_additive) %>%
  group_by(DKD_IND_additive) %>%
  dplyr::summarise(age_mean=mean(NVAL_NUM,na.rm=T),
                   age_sd=sd(NVAL_NUM,na.rm=T)) %>%
  ungroup %>% View

#other categorical
demo2<-demo %>%
  filter(!grepl("RELIGION",CONCEPT_CD)&NVAL_NUM>0) %>%
  dplyr::mutate(CONCEPT_TYPE=gsub("_.*","",CONCEPT_CD)) %>%
  group_by(DKD_IND_additive,CONCEPT_TYPE) %>%
  dplyr::mutate(enc_all=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  group_by(DKD_IND_additive,CONCEPT_TYPE,enc_all,CONCEPT_CD) %>%
  dplyr::summarize(enc_at=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  mutate(enc_prop=round(enc_at/enc_all,3))

demo2 %>%
  filter(!grepl("\\|",CONCEPT_CD)) %>%
  dplyr::select(CONCEPT_TYPE,CONCEPT_CD,enc_at,DKD_IND_additive) %>%
  spread(DKD_IND_additive,enc_at,fill=0) %>%
  left_join(demo2 %>%
              dplyr::select(CONCEPT_TYPE,CONCEPT_CD,enc_prop,DKD_IND_additive) %>%
              spread(DKD_IND_additive,enc_prop,fill=0),
            by=c("CONCEPT_TYPE","CONCEPT_CD")) %>%
  dplyr::rename("enc_cnt_0"=`0.x`,
                "enc_cnt_1"=`1.x`,
                "enc_prop_0"=`0.y`,
                "enc_prop_1"=`1.y`) %>%
  dplyr::select(CONCEPT_TYPE,CONCEPT_CD,
                enc_cnt_0,enc_prop_0,
                enc_cnt_1,enc_prop_1) %>%
  arrange(CONCEPT_TYPE,desc(enc_cnt_0)) %>%
  View


#dkd drift
dkd_drift<-pat_tbl %>%
  dplyr::select(PATIENT_NUM, DKD_IND_additive,YR_SINCE,tr_ts) %>%
  group_by(PATIENT_NUM,YR_SINCE) %>%
  top_n(n=1L,wt=DKD_IND_additive) %>%
  ungroup %>% unique %>%
  spread(YR_SINCE,DKD_IND_additive) %>%
  gather(YR_SINCE,DKD_IND_additive,-PATIENT_NUM,-tr_ts) %>%
  group_by(PATIENT_NUM) %>%
  mutate(DKD_IND=DKD_IND_additive) %>%
  fill(DKD_IND_additive,.direction="up") %>%
  ungroup %>%
  filter(!is.na(DKD_IND_additive)) %>%
  replace_na(list(DKD_IND="NA")) %>%
  dplyr::select(PATIENT_NUM,YR_SINCE,DKD_IND,tr_ts)

dkd_drift_demo<-dkd_drift %>%
  filter(DKD_IND != "NA") %>%
  left_join(demo %>% 
              filter(!grepl("\\|",CONCEPT_CD)) %>%
              dplyr::select(PATIENT_NUM,CONCEPT_CD,NVAL_NUM) %>%
              unique,
            by="PATIENT_NUM")


dkd_drift %>%
  filter(as.numeric(YR_SINCE)<=5) %>%
  group_by(DKD_IND,YR_SINCE,tr_ts) %>%
  dplyr::summarise(pat_cnt=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  spread(DKD_IND,pat_cnt,fill=0) %>%
  mutate(DKD_rt=round(`1`/(`0`+`1`),2)) %>%
  left_join(dkd_drift_demo %>%
              filter(CONCEPT_CD %in% c("AGE_at_DM","RACE_white","SEX_MALE","RELIGION_none")) %>%
              group_by(YR_SINCE,CONCEPT_CD,tr_ts) %>%
              dplyr::summarize(avg=mean(NVAL_NUM,na.rm=T)) %>%
              ungroup %>%
              mutate(CONCEPT_CD=paste0(CONCEPT_CD,"_avg")) %>%
              spread(CONCEPT_CD,avg),
            by=c("YR_SINCE","tr_ts")) %>%
  left_join(dkd_drift_demo %>%
              filter(CONCEPT_CD %in% c("AGE_at_DM","RACE_white","SEX_MALE","RELIGION_none")) %>%
              group_by(YR_SINCE,CONCEPT_CD,tr_ts) %>%
              dplyr::summarize(sd=sd(NVAL_NUM,na.rm=T)) %>%
              ungroup %>%
              mutate(CONCEPT_CD=paste0(CONCEPT_CD,"_sd")) %>%
              spread(CONCEPT_CD,sd),
            by=c("YR_SINCE","tr_ts")) %>%
  dplyr::mutate(eligb=`0`+`1`) %>%
  View



#heatmap for facts/patient over years
heatmap<-fact_stack %>%
  filter(yr_from_dm<=5) %>%
  mutate(yr_from_dm=case_when(yr_from_dm<=-5~-5,
                              TRUE~yr_from_dm)) %>%
  group_by(VARIABLE_CATEG,yr_from_dm) %>%
  dplyr::summarise(fact_cnt=n(),
                   pat_cnt=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  mutate(fact_per_pat=round(fact_cnt/pat_cnt))
saveRDS(heatmap,file="./output/heatmap.rda")

ggplot(heatmap %>% filter(!VARIABLE_CATEG %in% c("engineered","demographics")),
       aes(x=yr_from_dm,y=VARIABLE_CATEG))+
  geom_tile(aes(fill=fact_per_pat),color="white")+
  scale_fill_gradient(low="white",high="steelblue",name="# of facts/patient")+
  scale_x_continuous(breaks=seq(-5,5),labels = seq(-5,5))+
  labs(x="# of Years from DM onset",y="Data Type")+
  theme_classic()


#pick out one random DKD patients
#--patient_num=70
ex<-list(pat_tbl=pat_tbl %>% 
           filter(PATIENT_NUM==70) %>%
           dplyr::select(-PATIENT_NUM),
         fact_tbl=fact_stack %>%
           filter(PATIENT_NUM==70) %>%
           dplyr::select(-PATIENT_NUM))

saveRDS(ex,file="./output/example.rda")

ex<-readRDS("./output/example.rda")

ex_fact<-ex$fact_tbl %>%
  dplyr::mutate(day_from_dm=ifelse(VARIABLE_CATEG=="demographics",0,day_from_dm)) %>%
  dplyr::select(CONCEPT_CD,VARIABLE_CATEG,day_from_dm) %>%
  mutate(color=3) %>%
  bind_rows(ex$pat_tbl %>% 
              dplyr::select(DKD_IND_additive,DAY_SINCE) %>%
              gather(CONCEPT_CD,NVAL_NUM,-DAY_SINCE) %>%
              dplyr::rename(day_from_dm=DAY_SINCE) %>%
              mutate(VARIABLE_CATEG="DKD",color=1-NVAL_NUM) %>%
              dplyr::select(CONCEPT_CD,VARIABLE_CATEG,day_from_dm,color)) %>%
  bind_rows(ex$pat_tbl %>% 
              dplyr::select(SERUM_CREATININE,ACR,URINE_MICROALBUMIN,DAY_SINCE) %>%
              gather(CONCEPT_CD,NVAL_NUM,-DAY_SINCE) %>%
              dplyr::rename(day_from_dm=DAY_SINCE) %>%
              mutate(VARIABLE_CATEG="labs") %>%
              mutate(color=3) %>%
              dplyr::select(CONCEPT_CD,VARIABLE_CATEG,day_from_dm,color)) %>%
  bind_rows(ex$pat_tbl %>% 
              dplyr::select(ENC_TYPE,DAY_SINCE) %>%
              bind_rows(ex$fact_tbl %>%
                          dplyr::select(day_from_dm) %>%
                          dplyr::rename(DAY_SINCE=day_from_dm) %>%
                          mutate(ENC_TYPE="OA")) %>%
              gather(CONCEPT_CD,NVAL_NUM,-DAY_SINCE) %>%
              dplyr::rename(day_from_dm=DAY_SINCE) %>%
              mutate(VARIABLE_CATEG="encounters") %>%
              mutate(color=3) %>%
              dplyr::select(CONCEPT_CD,VARIABLE_CATEG,day_from_dm,color)) %>%
  filter(VARIABLE_CATEG!="engineered") %>%
  mutate(VARIABLE_CATEG=toupper(VARIABLE_CATEG))

ex_fact$VARIABLE_CATEG<-factor(ex_fact$VARIABLE_CATEG,
                               levels=c("DKD",
                                        "VISIT_DETAILS",
                                        "ORDERS",
                                        "MEDICATIONS",
                                        "LABS",
                                        "HISTORY",
                                        "DIAGNOSES",
                                        "DEMOGRAPHICS",
                                        "ALLERGY",
                                        "ALERTS",
                                        "ENCOUNTERS"))

ggplot(ex_fact %>% filter(day_from_dm>=-365),
       aes(x=day_from_dm,y=VARIABLE_CATEG,color=as.factor(color)))+
  geom_point(alpha=0.3,size=10)+
  scale_shape(solid=F)+
  scale_color_discrete(guide=F)+
  scale_x_continuous("Days since DM onset",
                     breaks=c(0:6)*360,
                     labels=paste0("+",0:6,"years"))+
  geom_vline(xintercept=c(0:6)*360,
             linetype=2, 
             color="black",
             size=1)+
  labs(y="Data Type")+
  theme_hc()

