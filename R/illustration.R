## illustrations ##

rm(list=ls())

source("./R/util.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
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



