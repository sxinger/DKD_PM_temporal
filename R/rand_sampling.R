#### random sample ####

rm(list=ls()); gc()

source("./R/util.R")
require_libraries(c("dplyr",
                    "tidyr",
                    "magrittr",
                    "plyr"
))


## Load targets
pat_tbl<-read.csv("./data/dkd_patients_full.csv",header=T,stringsAsFactors=F) %>%
  semi_join(read.csv("./data/dkd_patients.csv",header=T,stringsAsFactors=F),
            by="PATIENT_NUM") %>%
  dplyr::mutate(ENC_START = as.Date(ENC_START,format="%d-%b-%y"),
                ENC_END = as.Date(ENC_END,format="%d-%b-%y"),
                DM_ONSET_DATE = as.Date(DM_ONSET_DATE,format="%d-%b-%y"),
                LAB_DATE = as.Date(LAB_DATE,format="%d-%b-%y"),
                END_DATE = as.Date(END_DATE,format="%d-%b-%y")) %>%
  group_by(PATIENT_NUM) %>% arrange(END_DATE) %>%
  dplyr::mutate(DKD_IND_additive = cumsum(DKD_IND)) %>% ungroup %>%
  dplyr::mutate(DM_YR = as.numeric(format(DM_ONSET_DATE,"%Y")),
                END_YR =  as.numeric(format(END_DATE,"%Y"))) %>%
  dplyr::mutate(DAY_SINCE = difftime(END_DATE,DM_ONSET_DATE,units="days")) %>%
  dplyr::mutate(DAY_SINCE = as.numeric(DAY_SINCE)) %>%
  group_by(PATIENT_NUM) %>%
  dplyr::mutate(DAY_SINCE_lag = lag(DAY_SINCE,n=1L,default = 0)) %>%
  dplyr::mutate(DAY_SINCE_delta = DAY_SINCE - DAY_SINCE_lag) %>%
  dplyr::mutate(MTH1_SINCE = floor(DAY_SINCE/30),
                MTH3_SINCE = floor(DAY_SINCE/90),
                MTH6_SINCE = floor(DAY_SINCE/182.5),
                YR_SINCE = floor(DAY_SINCE/365.25)) %>%
  arrange(PATIENT_NUM, ENCOUNTER_NUM, ENC_START) %>%
  filter(DKD_IND_additive <=1 & YR_SINCE >=0) %>%
  ungroup

## Simple Random Sampling
pn_fold<-pat_tbl %>% 
  dplyr::select(PATIENT_NUM) %>% unique %>%
  dplyr::mutate(part_idx=sample(1:5,n(),replace=T)) 

pat_tbl %<>% 
  left_join(pn_fold,by="PATIENT_NUM")

saveRDS(pat_tbl,file="./data2/pat_episode2.rda")

