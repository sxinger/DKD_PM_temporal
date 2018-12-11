## illustrations ##

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


load("./data/DKD_heron_facts_prep.Rdata")
load("./data2/pat_episode2.Rdata")

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
                   cd_cnt=length(unique(CONCEPT_CD))) %>%
  ungroup

#heatmap for facts/patient over years
heatmap<-fact_stack %>%
  mutate(yr_from_dm=ifelse(yr_from_dm<=-5,-5,yr_from_dm)) %>%
  group_by(VARIABLE_CATEG,yr_from_dm) %>%
  dplyr::summarise(fact_cnt=n(),
                   pat_cnt=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  mutate(fact_per_pat=round(fact_cnt/pat_cnt))
save(heatmap,file="./output/heatmap.rda")

df<-fact_stack %>%
  group_by(VARIABLE_CATEG) %>%
  dplyr::summarise(fact_cnt=n(),
                   pat_cnt=length(unique(PATIENT_NUM))) %>%
  ungroup %>%
  mutate(fact_per_pat=round(fact_cnt/pat_cnt))


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

