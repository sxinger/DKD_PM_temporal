#### Feature Engineering ####

rm(list=ls()); gc()

setwd("~/proj_dkd/DKD_PM_wip")

source("./util.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
))

fact_stack<-readRDS("./data2/DKD_heron_facts_prep.rda")

## add historical distinct fact counts (distinct day,concept) update to date
add_FCNT<-fact_stack %>% 
  dplyr::select(PATIENT_NUM, CONCEPT_CD, day_from_dm) %>%
  dplyr::mutate(day_from_dm = pmax(0,day_from_dm)) %>% unique %>%
  group_by(PATIENT_NUM,day_from_dm) %>% 
  dplyr::mutate(fcnt_pday = n()) %>%
  ungroup %>% arrange(PATIENT_NUM,day_from_dm) %>%
  dplyr::select(PATIENT_NUM,day_from_dm,fcnt_pday) %>% unique %>%
  dplyr::mutate(VARIABLE_CATEG = "engineered",
                CONCEPT_CD = "fact_cnt",
                NVAL_NUM = cumsum(fcnt_pday)) %>%
  dplyr::mutate(NVAL_NUM_lag = lag(NVAL_NUM,n=1L,0)) %>%
  dplyr::mutate(NVAL_NUM2 = NVAL_NUM-NVAL_NUM_lag)


add_FCNT2<-fact_stack %>%
  dplyr::select(-VARIABLE_CATEG,-CONCEPT_CD,-NVAL_NUM) %>% unique %>%
  inner_join(add_FCNT %>% dplyr::mutate(CONCEPT_CD = "fact_cnt") %>%
             dplyr::select(PATIENT_NUM,day_from_dm,
                            VARIABLE_CATEG,CONCEPT_CD,NVAL_NUM),
            by=c("PATIENT_NUM","day_from_dm")) %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, CONCEPT_CD,
                START_DATE,DM_DATE,END_DATE,
                day_from_dm,day_to_end,yr_from_dm,NVAL_NUM)

#eyeball an example
add_FCNT2 %>% filter(PATIENT_NUM==70) %>% View


add_FCNT3<-fact_stack %>%
  dplyr::select(-VARIABLE_CATEG,-CONCEPT_CD,-NVAL_NUM) %>% unique %>%
  inner_join(add_FCNT %>% dplyr::mutate(CONCEPT_CD = "newfact_cnt_slast") %>%
              dplyr::select(PATIENT_NUM,day_from_dm,
                            VARIABLE_CATEG,CONCEPT_CD,NVAL_NUM2) %>%
              dplyr::rename(NVAL_NUM = NVAL_NUM2),
            by=c("PATIENT_NUM","day_from_dm")) %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, CONCEPT_CD, START_DATE,DM_DATE,END_DATE,
                day_from_dm,day_to_end,yr_from_dm,NVAL_NUM)

#eyeball an example
add_FCNT3 %>% filter(PATIENT_NUM==70) %>% View


## stack new features
fact_stack %<>%
  bind_rows(add_FCNT2)
rm(add_FCNT2); gc()

fact_stack %<>%
  bind_rows(add_FCNT3)
rm(add_FCNT3); gc()

## update data
saveRDS(fact_stack, file="./data2/DKD_heron_facts_prep.rda")

rm(list=ls())
gc()
