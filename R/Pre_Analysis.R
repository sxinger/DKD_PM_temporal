#### Pre-Analysis ####
rm(list=ls())
gc()

source("./R/util.R")
require_libraries(c( "Matrix"
                     ,"pROC"
                     ,"xgboost"
                     ,"dplyr"
                     ,"tidyr"
                     ,"magrittr"
                   ))

pat_tbl<-readRDS("./data2/pat_episode2.rda")
fact_stack<-readRDS("./data2/DKD_heron_facts_prep.rda")

#### eGFR update frequencies ####
eGFR_ud_freq<-pat_tbl %>%
  group_by(PATIENT_NUM) %>%
  dplyr::summarize(day_delta_mean = mean(DAY_SINCE_delta,na.rm=T),
                   day_delta_sd = sd(DAY_SINCE_delta,na.rm=T),
                   day_delta_median = median(DAY_SINCE_delta,na.rm=T),
                   day_delta_IQR = quantile(DAY_SINCE_delta,probs=0.75,na.rm=T)-quantile(DAY_SINCE_delta,probs=0.25,na.rm=T)) %>%
  ungroup %>%
  dplyr::summarize(overall_mean = mean(day_delta_mean,na.rm=T),
                   within_pat_sd = mean(day_delta_sd,na.rm=T),
                   acr_pat_sd = sd(day_delta_mean,na.rm=T),
                   overall_median = median(day_delta_median,na.rm=T),
                   within_pat_IQR = mean(day_delta_IQR,na.rm=T),
                   acr_pat_IQR = sd(day_delta_IQR,na.rm=T))


#### clinical fact intensity ####
fact_ud_freq<-fact_stack %>%
  dplyr::select(PATIENT_NUM,VARIABLE_CATEG,day_from_dm) %>% 
  dplyr::mutate(day_from_dm = pmax(0,day_from_dm)) %>%
  unique %>%
  group_by(PATIENT_NUM,VARIABLE_CATEG) %>%
  arrange(day_from_dm) %>%
  dplyr::mutate(day_from_dm_lag = lag(day_from_dm,n=1L)) %>%
  dplyr::mutate(day_from_dm_delta = day_from_dm - day_from_dm_lag) %>%
  dplyr::summarize(day_delta_mean = mean(day_from_dm_delta,na.rm=T),
                   day_delta_sd = sd(day_from_dm_delta,na.rm=T),
                   day_delta_median = median(day_from_dm_delta,na.rm=T),
                   day_delta_IQR = quantile(day_from_dm_delta,probs=0.75,na.rm=T)-quantile(day_from_dm_delta,probs=0.25,na.rm=T)) %>%
  ungroup %>%
  group_by(VARIABLE_CATEG) %>%
  dplyr::summarize(size=length(unique(PATIENT_NUM)),
                   overall_mean = mean(day_delta_mean,na.rm=T),
                   within_pat_sd = mean(day_delta_sd,na.rm=T),
                   acr_pat_sd = sd(day_delta_mean,na.rm=T),
                   overall_median = median(day_delta_median,na.rm=T),
                   within_pat_IQR = mean(day_delta_IQR,na.rm=T),
                   acr_pat_IQR = sd(day_delta_IQR,na.rm=T))
