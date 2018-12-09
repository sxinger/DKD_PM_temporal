#### Preprocessing ####

rm(list=ls()); gc()

setwd("~/proj_dkd/DKD_PM_wip")

source("./util.R")
require_libraries(c("Matrix",
                    "dplyr",
                    "tidyr",
                    "plyr",
                    "magrittr",
                    # "stringr",
                    "moments"
))

## patient drift 
pat_tbl<-read.csv("./data/dkd_patients_full.csv",header=T,stringsAsFactors=F) %>%
  dplyr::mutate(ENC_START = as.Date(ENC_START,format="%d-%b-%y"),
                ENC_END = as.Date(ENC_END,format="%d-%b-%y"),
                DM_DATE = as.Date(DM_ONSET_DATE,format="%d-%b-%y"),
                LAB_DATE = as.Date(LAB_DATE,format="%d-%b-%y"),
                END_DATE = as.Date(END_DATE,format="%d-%b-%y")) %>%
  group_by(PATIENT_NUM) %>% arrange(END_DATE) %>%
  dplyr::mutate(DKD_IND_additive = cumsum(DKD_IND)) %>%
  dplyr::mutate(DM_YR = as.numeric(format(DM_DATE,"%Y")),
                END_YR =  as.numeric(format(END_DATE,"%Y"))) %>%
  dplyr::mutate(DAY_SINCE = difftime(END_DATE,DM_DATE,units="days")) %>%
  dplyr::mutate(YR_SINCE = round(DAY_SINCE/365.25)) %>%
  arrange(PATIENT_NUM, END_DATE) %>%
  filter(YR_SINCE >= 0) %>% filter(DKD_IND_additive <=1) %>%
  ungroup

# average eGFR over onset years
eGFR_over_yr<-pat_tbl %>%
  dplyr::filter(DM_YR > 2006) %>%
  select(PATIENT_NUM,DM_ONSET_DATE,DM_YR,END_DATE,END_YR,EGFR_MDRD,YR_SINCE) %>% unique %>%
  group_by(DM_YR,YR_SINCE) %>%
  dplyr::summarize(pat_cnt=length(unique(PATIENT_NUM)),
                   eGFR_expos=sum(!is.na(EGFR_MDRD)),
                   min_eGFR=min(EGFR_MDRD,na.rm=T),
                   q1_eGFR=quantile(EGFR_MDRD,probs=0.25,na.rm=T),
                   mean_eGFR=mean(EGFR_MDRD,na.rm=T),
                   median_eGFR=median(EGFR_MDRD,na.rm=T),
                   q3_eGFR=quantile(EGFR_MDRD,probs=0.75,na.rm=T),
                   max_eGFR=max(EGFR_MDRD,na.rm=T),
                   eGFR_less60=mean((EGFR_MDRD<60),na.rm=T))

eGFR_over_yr<-pat_tbl %>%
  dplyr::filter(DM_YR > 2006) %>%
  select(PATIENT_NUM,DM_ONSET_DATE,DM_YR,END_DATE,END_YR,EGFR_IN_EHR,YR_SINCE) %>% unique %>%
  group_by(DM_YR,YR_SINCE) %>%
  dplyr::summarize(pat_cnt=length(unique(PATIENT_NUM)),
                   eGFR_expos=sum(!is.na(EGFR_IN_EHR)),
                   min_eGFR=min(EGFR_IN_EHR,na.rm=T),
                   q1_eGFR=quantile(EGFR_IN_EHR,probs=0.25,na.rm=T),
                   mean_eGFR=mean(EGFR_IN_EHR,na.rm=T),
                   median_eGFR=median(EGFR_IN_EHR,na.rm=T),
                   q3_eGFR=quantile(EGFR_IN_EHR,probs=0.75,na.rm=T),
                   max_eGFR=max(EGFR_IN_EHR,na.rm=T),
                   eGFR_less60=mean((EGFR_IN_EHR<60),na.rm=T)
  )

Scr_over_yr<-pat_tbl %>%
  dplyr::filter(DM_YR > 2006) %>%
  select(PATIENT_NUM,DM_ONSET_DATE,DM_YR,END_DATE,END_YR,SERUM_CREATININE,YR_SINCE) %>% unique %>%
  group_by(DM_YR,YR_SINCE) %>%
  dplyr::summarize(pat_cnt=length(unique(PATIENT_NUM)),
                   Scr_expos=sum(!is.na(SERUM_CREATININE)),
                   min_Scr=min(SERUM_CREATININE,na.rm=T),
                   q1_Scr=quantile(SERUM_CREATININE,probs=0.25,na.rm=T),
                   mean_Scr=mean(SERUM_CREATININE,na.rm=T),
                   median_Scr=median(SERUM_CREATININE,na.rm=T),
                   q3_Scr=quantile(SERUM_CREATININE,probs=0.75,na.rm=T),
                   max_Scr=max(SERUM_CREATININE,na.rm=T)
  )

# Load patient table
pat_tbl<-read.csv("./data/dkd_patients.csv",header=T,stringsAsFactors=F) %>%
  dplyr::mutate(BIRTH_DATE = as.Date(BIRTH_DATE,"%d-%m-%Y %H:%M:%S"),
                DEATH_DATE = as.Date(DEATH_DATE,"%d-%m-%Y %H:%M:%S"),
                DM_DATE = as.Date(DM_DATE,"%d-%m-%Y %H:%M:%S"),
                END_DATE = as.Date(END_DATE,"%d-%m-%Y %H:%M:%S")) %>%
  dplyr::mutate(DAY_BT_BIRTH_DM = as.numeric(difftime(DM_DATE,BIRTH_DATE,units="days")),
                DAY_SINCE = as.numeric(difftime(END_DATE,DM_DATE,units="days")),
                DAY_TO = as.numeric(difftime(END_DATE,BIRTH_DATE,units="days"))) %>%
  dplyr::mutate(DM_YR = as.numeric(format(DM_DATE,"%Y")),
                END_YR =  as.numeric(format(END_DATE,"%Y"))) %>%
  dplyr::mutate(YR_SINCE = round(DAY_SINCE/365.25),
                YR_TO = round(DAY_TO/365.25),
                AGE_at_DM = round(DAY_BT_BIRTH_DM/365.25)) %>%
  dplyr::select(PATIENT_NUM, BIRTH_DATE, DM_DATE, END_DATE, DEATH_DATE,
                DM_YR, END_YR, DAY_SINCE, DAY_TO, YR_SINCE, YR_TO,
                DKD_IND, AGE_at_DM, SEX_MALE, RACE, RELIGION)

save(pat_tbl, file="./data/pat_tbl.Rdata")


## Load in original data sets and only track value change on yearly basis
data_type<-c("allergy","demographics",
             "history","uhc","visit_details",
             "alerts","diagnoses","procedures","orders","medications",
             "labs")

fact_stack<-c()
for(i in 1:length(data_type)){
  dat_i<-read.csv(paste0("./data/dkd_",data_type[i],".csv"),header=T,stringsAsFactors=F) %>%
    dplyr::select(PATIENT_NUM, CONCEPT_CD, NVAL_NUM, START_DATE, RECENCY, C_NAME, C_VISUAL_PATH) %>%
    dplyr::mutate(START_DATE=as.Date(START_DATE,format="%d-%b-%y"),
                  VARIABLE_CATEG=data_type[i])
  fact_stack %<>% bind_rows(dat_i)
}
#stack additional demographic facts collected from patient_dimension
fact_stack %<>%
  bind_rows(pat_tbl %>% 
              dplyr::select(PATIENT_NUM,BIRTH_DATE,AGE_at_DM,SEX_MALE,RACE,RELIGION) %>%
              mutate(RACE=paste0("RACE_",gsub(" ","_",RACE)),
                     RELIGION=paste0("RELIGION_",gsub(" ","_",RELIGION))) %>%
              mutate(race_ref=1, relig_ref=1) %>%
              spread(RACE, race_ref, fill=0) %>%
              spread(RELIGION, relig_ref, fill=0) %>%
              gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM,-BIRTH_DATE) %>%
              mutate(START_DATE = as.Date(BIRTH_DATE,"%Y-%m-%d"),
                     RECENCY = 1,
                     C_NAME = CONCEPT_CD,
                     C_VISUAL_PATH = "patient_dimension",
                     VARIABLE_CATEG = "demographics") %>%
              select(-BIRTH_DATE))
# save raw data
save(fact_stack,file="./data/fact_stack_episode.Rdata")


## Pre-Process
load("./data/fact_stack_episode.Rdata")
unique(fact_stack$VARIABLE_CATEG)
# [1] "allergy"       "demographics"  "history"       "uhc"           "visit_details" "alerts"        "diagnoses"     "procedures"    "orders"        "medications"  
# [11] "labs"

## Preprocess 1 --- pre-select data subtypes
# - selected concepts from demographic
cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  dplyr::filter(!(VARIABLE_CATEG == "demographics" & !grepl("(ETHNICITY)|(GEO)|(SEX)|(RACE)",CONCEPT_CD)))
cnt_after<-nrow(fact_stack)
#report 
cat("Removed",cnt_before - cnt_after,"facts.\n")
# Removed 15,421 facts.


##Preprocess 2 --- Impute NULL
nval_summary<-fact_stack %>%
  group_by(VARIABLE_CATEG,CONCEPT_CD) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   distinct_val_narm=length(unique(NVAL_NUM))-ifelse(sum(is.na(NVAL_NUM))>0,1,0),
                   na_prop=round(sum(is.na(NVAL_NUM))/n(),4),
                   neg_prop=round(sum((NVAL_NUM<0))/n(),4),
                   zero_prop=round(sum((NVAL_NUM==0))/n(),4),
                   min_val=ifelse(sum(is.na(NVAL_NUM))/n()==1,NA,min(NVAL_NUM,na.rm=T)),
                   max_val=ifelse(sum(is.na(NVAL_NUM))/n()==1,NA,max(NVAL_NUM,na.rm=T)))

nval_curation<-list(neg_val=nval_summary %>% filter(neg_prop>0),
                    num_val_w_na=nval_summary %>% filter(distinct_val_narm>0 & na_prop>0))
save(nval_curation,file="./data/nval_num_abnormality.Rdata")

nval_summary %<>% ungroup %>%
  select(CONCEPT_CD,distinct_val_narm)

na_cnt_before<-nrow(fact_stack %>% filter(is.na(NVAL_NUM)))
fact_stack %<>%
  left_join(nval_summary,by="CONCEPT_CD") %>%
  mutate(NVAL_NUM=ifelse(distinct_val_narm==0 & is.na(NVAL_NUM),1,
                         ifelse(is.na(NVAL_NUM),0,NVAL_NUM)))
na_cnt_after<-nrow(fact_stack %>% filter(is.na(NVAL_NUM)))
#report 
cat("Impute",na_cnt_before - na_cnt_after,"values for nval_num \n")
# Impute 18,377,361 values for nval_num (if binary, impute 1; if numeric, impute 0)

#sanity checks
# summary(fact_stack %>% filter(CONCEPT_CD=="KUMC|PACK_PER_DAY") %>% 
#         select(NVAL_NUM) %>% unlist)
# 
# summary(fact_stack %>% filter(CONCEPT_CD=="KUH|COMPONENT_ID:3017#%@Labs|Aggregate:Median") %>% 
#         select(NVAL_NUM) %>% unlist)
# 
# summary(fact_stack %>% filter(CONCEPT_CD=="KUH|BPA_LOCATOR_ID:1029@ALT_STATUS_C:0") %>% 
#         select(NVAL_NUM) %>% unlist)

# View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
#        dplyr::summarize(fact_cnt=n(),
#                         pat_cnt=length(unique(PATIENT_NUM)),
#                         cd_cnt=length(unique(CONCEPT_CD))))

## Preprocess 3 --- Correct possible typos
# - temperature conversion (C to F)
fact_stack %<>%
  mutate(NVAL_NUM = ifelse(CONCEPT_CD=="KUH|PAT_ENC:TEMPERATURE" & NVAL_NUM<47,NVAL_NUM*1.8+32,NVAL_NUM)) 
gc()


## Preprocess 4 --- Identify and remove possible outliers
# using prior knowledge
outliers<-read.csv("./data/variable_common_range.csv",stringsAsFactors = F)

cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  left_join(outliers, by="CONCEPT_CD") %>%
  dplyr::filter(is.na(upper)|(NVAL_NUM <= upper & NVAL_NUM > lower)) %>%
  dplyr::select(-lower,-upper)
cnt_after<-nrow(fact_stack)
#report 
cat("Identified and removed",cnt_before - cnt_after,"facts.\n")
# Identified and removed 17 facts.
gc()


## Preprocess 5 -- take out medication and diagnoses facts
heron_terms<-read.csv("./data/heron_terms.csv",header=T,stringsAsFactors=F)
unique(heron_terms$VARIABLE_CATEG)
# [1] "CF"                 "PROCEDURES"         "NAACCR"             "DEMOGRAPHICS"       "VISIT DETAILS"      "DIAGNOSES"          "ALLERGY"           
# [8] "ALERTS"             "LABORATORY TESTS"   "LABTESTS"           "CARDIOLABTESTS"     "FLOWSHEET"          "MEDICATIONS"        "PROCORDERS"        
# [15] "HISTORY"            "SPECIMENS"          "MICRONEGATIVE"      "MICROPOSITIVE"      "MICROBIOLOGY"       "REPORTS"            "RESEARCHENROLLMENT"
# [22] "ABRIDGED"           "NCDR"               "NTDS"               "UHC"           

med_stack<-fact_stack %>%
  dplyr::filter(VARIABLE_CATEG == "medications") %>%
  mutate(link_concept=gsub("[#@].*$","",CONCEPT_CD)) %>%
  left_join(heron_terms %>% filter(VARIABLE_CATEG=="MEDICATIONS") %>%
              dplyr::select(C_BASECODE, C_NAME, C_VISUAL_PATH) %>% unique,
            by=c("link_concept"="C_BASECODE")) %>%
  select(-link_concept)

dx_stack<-fact_stack %>%
  dplyr::filter(VARIABLE_CATEG == "diagnoses") %>%
  mutate(link_concept=gsub("[#@].*$","",CONCEPT_CD)) %>%
  left_join(heron_terms %>% filter(VARIABLE_CATEG=="DIAGNOSES") %>%
              dplyr::select(C_BASECODE, C_NAME, C_VISUAL_PATH) %>% unique,
            by=c("link_concept"="C_BASECODE")) %>%
  select(-link_concept)

other_fact_stack<-fact_stack %>%
  dplyr::filter(!VARIABLE_CATEG %in% c("medications","diagnoses"))


## Preprocess 6 -- feature engineering
#medication - along the HERON-ontology tree
#-- raw value
med_stack %<>%
  mutate(C_VISUAL_PATH=gsub("\\\\i2b2\\\\","",C_VISUAL_PATH)) %>%
  mutate(NLEVEL=stringr::str_count(C_VISUAL_PATH,"\\\\"),
         VA_NLEVEL=stringr::str_count(C_VISUAL_PATH,"\\[")) %>%
  mutate(RX_MODIFIER=gsub("^.*[@]","",CONCEPT_CD),
         CONCEPT_LEVEL=6)

med_stack2<-med_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                RX_MODIFIER, CONCEPT_LEVEL)

#-- concept level
med_stack %<>%
  mutate(CONCEPT_CD=gsub("[#@].*$","",CONCEPT_CD)) %>%
  mutate(NVAL_NUM=1, CONCEPT_LEVEL=5)

med_stack_new<-med_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                RX_MODIFIER, CONCEPT_LEVEL)
med_stack2 %<>%
  bind_rows(med_stack_new)

#--SCDF or SCBF
med_stack_new<-med_stack %>%
  mutate(C_VISUAL_PATH = ifelse(NLEVEL==(VA_NLEVEL+1),
                                C_VISUAL_PATH,
                                ifelse(VA_NLEVEL==0,
                                       stringr::str_extract(C_VISUAL_PATH,
                                                            paste0("^([^\\\\]*\\\\){",pmin(3,NLEVEL),"}")),
                                       stringr::str_extract(C_VISUAL_PATH,
                                                            paste0("^([^\\\\]*\\\\){",(VA_NLEVEL+2),"}"))))) %>%
  mutate(C_NAME=gsub("[^\\\\]*\\\\","",gsub("\\\\$","",C_VISUAL_PATH))) %>%
  mutate(CONCEPT_CD=C_NAME, NVAL_NUM=1, CONCEPT_LEVEL=4) %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                RX_MODIFIER, CONCEPT_LEVEL)
med_stack2 %<>%
  bind_rows(med_stack_new)

# # #eyeball random dample
# View(med_stack_new[sample(seq_len(nrow(med_stack_new)),50),])

#--VA class levels
va_max<-max(med_stack$VA_NLEVEL)
for(lev in 1:va_max){
  med_stack_new<-med_stack %>%
    mutate(C_VISUAL_PATH=ifelse(VA_NLEVEL==0, 
                                stringr::str_extract(C_VISUAL_PATH,
                                                     paste0("^([^\\\\]*\\\\){2}")),
                                stringr::str_extract(C_VISUAL_PATH,
                                                     paste0("^([^\\\\]*\\\\){",
                                                            pmin((VA_NLEVEL+1),(lev+1)),
                                                            "}")))) %>%
    # mutate(RX_CLASS_UD=stringr::str_extract(C_VISUAL_PATH,
    #                                         paste0("^([^\\\\]*\\\\){",
    #                                                NLEVEL-pmin(lev,(VA_NLEVEL-1))+1,
    #                                                "}"))) %>%
    mutate(C_NAME=gsub("[^\\\\]*\\\\","",gsub("\\\\$","",C_VISUAL_PATH))) %>%
    mutate(CONCEPT_CD=coalesce(gsub("\\[|\\]","",stringr::str_extract(C_NAME,"\\[.*?\\]")),C_NAME)) %>%
    mutate(NVAL_NUM=1, CONCEPT_LEVEL=lev) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                  CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                  RX_MODIFIER, CONCEPT_LEVEL)
  
  med_stack2 %<>%
    bind_rows(med_stack_new)
}

#eyeball random dample
# View(med_stack2[sample(seq_len(nrow(med_stack2)),100),])
# med_stack2 %>% filter(PATIENT_NUM==388225 & START_DATE=="2016-02-15" & RX_MODIFIER=="Surescripts:Amount") %>% unique %>% View

#diagnoses - along the HERON-ontology tree
#-- at raw level
dx_stack %<>%
  mutate(C_VISUAL_PATH=stringr::str_replace(C_VISUAL_PATH,"\\\\i2b2\\\\","")) %>%
  dplyr::filter(!grepl("^(Diagnoses\\\\A)+", C_VISUAL_PATH)) %>% #remove
  mutate(C_VISUAL_PATH=paste0(stringr::str_replace(C_VISUAL_PATH,"Diagnoses\\\\",""),"\\")) %>%
  mutate(DX_type = stringr::str_replace(stringr::str_extract(C_VISUAL_PATH,"^([^\\\\]*\\\\){1}"),"\\\\","")) %>%
  mutate(NLEVEL=stringr::str_count(C_VISUAL_PATH,"\\\\")) %>%
  mutate(DX_MODIFIER=gsub("^.*[@]","",CONCEPT_CD),
         CONCEPT_LEVEL=7)

dx_stack2<-dx_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                DX_MODIFIER, CONCEPT_LEVEL)

#-- at concept level
dx_stack %<>% 
  mutate(CONCEPT_CD=gsub("[@].*$","",CONCEPT_CD),
         CONCEPT_LEVEL=6) %>%
  mutate(DECIMAL=ifelse(grepl("^(Other Diagnoses Concepts)+",C_VISUAL_PATH), 1,
                        ifelse(grepl("DX_ID",CONCEPT_CD),
                               pmax(0,
                                    nchar(stringr::str_replace(paste0(stringr::str_replace(stringr::str_replace(C_VISUAL_PATH,
                                                                                     paste0("^([^\\\\]*\\\\){",(NLEVEL-2),"}"),""),
                                                                         "\\s.*",""),
                                                             "."),
                                          "[^\\.]*(\\.)",""))-1) + 1,
                               pmax(0,
                                    nchar(stringr::str_replace(paste0(CONCEPT_CD,"."),"[^\\.]*(\\.)","")) - 1))))

dx_stack_new<-dx_stack %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                DX_MODIFIER, CONCEPT_LEVEL)
dx_stack2 %<>%
  bind_rows(dx_stack_new)


#-- below DX interger level
# dx_stack %>% dplyr::filter(DECIMAL>3 & grepl("DX_ID",CONCEPT_CD)) 

for(lev in 4:5){
  dx_stack_new<-dx_stack %>%
    mutate(C_VISUAL_PATH=ifelse(grepl("DX_ID",CONCEPT_CD),
                                stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",
                                                                          pmin((NLEVEL-1),((NLEVEL-DECIMAL)+(lev-3))),"}")),
                                stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",
                                                                          pmin(NLEVEL,((NLEVEL-DECIMAL)+(lev-3))),"}")))) %>%
    # mutate(C_VISUAL_PATH = stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",NLEVEL-pmin(lev,(NLEVEL-2))+1,"}"))) %>%
    mutate(C_NAME=stringr::str_replace_all(gsub("\\\\$","",C_VISUAL_PATH),"[^\\\\]*\\\\","")) %>%
    mutate(CONCEPT_CD=paste0(DX_type,":",
                             gsub("\\(|\\)","",stringr::str_replace(C_NAME,"\\s.*","")))) %>%
    mutate(CONCEPT_LEVEL=lev) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                  CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                  DX_MODIFIER, CONCEPT_LEVEL)
  
  dx_stack2 %<>%
    bind_rows(dx_stack_new)
}


#--at DX integer level
dx_stack_new<-dx_stack %>%
  mutate(C_VISUAL_PATH=stringr::str_extract(C_VISUAL_PATH,
                                            paste0("^([^\\\\]*\\\\){",
                                                   (NLEVEL-DECIMAL),"}"))) %>%
  # mutate(C_VISUAL_PATH = stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",(NLEVEL-DECIMAL),"}"))) %>%
  mutate(C_NAME=stringr::str_replace_all(gsub("\\\\$","",C_VISUAL_PATH),"[^\\\\]*\\\\","")) %>%
  mutate(CONCEPT_CD=paste0(DX_type,":",
                           gsub("\\(|\\)","",stringr::str_replace(C_NAME,"\\s.*","")))) %>%
  mutate(CONCEPT_LEVEL=3) %>%
  dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                DX_MODIFIER, CONCEPT_LEVEL)

dx_stack2 %<>%
  bind_rows(dx_stack_new)


#--above DX interger level
# dx_stack %>% dplyr::filter(NLEVEL-DECIMAL>4) %>% head

for(lev in 1:2){
  dx_stack_new<-dx_stack %>%
    mutate(C_VISUAL_PATH=stringr::str_extract(C_VISUAL_PATH,
                                     paste0("^([^\\\\]*\\\\){",pmin((NLEVEL-DECIMAL-1),(lev+1)),"}"))) %>%
    # mutate(C_VISUAL_PATH = stringr::str_extract(C_VISUAL_PATH,paste0("^([^\\\\]*\\\\){",NLEVEL-pmin(lev,(NLEVEL-2))+1,"}"))) %>%
    mutate(C_NAME=stringr::str_replace(gsub("\\\\$","",C_VISUAL_PATH),"[^\\\\]*\\\\","")) %>%
    mutate(CONCEPT_CD=paste0(DX_type,":",
                             gsub("\\(|\\)","",stringr::str_replace(C_NAME,"\\s.*","")))) %>%
    mutate(CONCEPT_LEVEL=lev) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_CATEG, C_VISUAL_PATH, 
                  CONCEPT_CD, C_NAME, NVAL_NUM, START_DATE, 
                  DX_MODIFIER, CONCEPT_LEVEL)
  
  dx_stack2 %<>%
    bind_rows(dx_stack_new)
}

#eyeball random dample
# View(dx_stack2[sample(seq_len(nrow(dx_stack2)),100),])
# dx_stack2 %>% filter(PATIENT_NUM==723795 & START_DATE=="2008-10-17" & DX_MODIFIER=="DiagObs:PAT_ENC_DX") %>% unique %>% View

#save tables
save(other_fact_stack, file="./data/DKD_heron_other_facts_prep.Rdata")
save(med_stack2, file="./data/DKD_heron_med_facts_prep.Rdata")
save(dx_stack2, file="./data/DKD_heron_dx_facts_prep.Rdata")


#### Inclusion of Medicaion and Diagnoses features ####
#if breaks:
load("./data/DKD_heron_other_facts_prep.Rdata")
load("./data/DKD_heron_med_facts_prep.Rdata")
load("./data/DKD_heron_dx_facts_prep.Rdata")

#update demographics
demo_from_obs<-other_fact_stack %>% 
  dplyr::filter((VARIABLE_CATEG == "demographics" & !grepl("((ETHNICITY)|(GEO))+",CONCEPT_CD)))

demo_from_pat<-pat_tbl %>% 
  dplyr::select(PATIENT_NUM,BIRTH_DATE,AGE_at_DM,SEX_MALE,RACE,RELIGION,
                DAY_SINCE, DAY_TO, YR_SINCE, YR_TO) %>%
  mutate(RACE=paste0("RACE_",gsub(" ","_",RACE)),
         RELIGION=paste0("RELIGION_",gsub(" ","_",RELIGION))) %>%
  mutate(race_ref=1, relig_ref=1) %>%
  spread(RACE, race_ref, fill=0) %>%
  spread(RELIGION, relig_ref, fill=0) %>%
  gather(CONCEPT_CD,NVAL_NUM,-PATIENT_NUM,-BIRTH_DATE,
         -DAY_SINCE, -DAY_TO, -YR_SINCE, -YR_TO) %>%
  mutate(START_DATE = as.Date(BIRTH_DATE,"%Y-%m-%d"),
         RECENCY = 1,
         C_NAME = CONCEPT_CD,
         C_VISUAL_PATH = "patient_dimension",
         VARIABLE_CATEG = "demographics",
         day_from_dm = DAY_SINCE,
         day_to_end = DAY_TO,
         yr_from_dm = YR_SINCE,
         yr_to_end = YR_TO) %>%
  dplyr::select(PATIENT_NUM, CONCEPT_CD, NVAL_NUM,
                START_DATE,RECENCY, VARIABLE_CATEG,
                day_from_dm, day_to_end,
                yr_from_dm, yr_to_end)

other_fact_stack %<>%
  filter(VARIABLE_CATEG != "demographics") %>%
  dplyr::select(-distinct_val_narm) %>%
  dplyr::mutate(day_from_dm = as.numeric(day_from_dm),
                day_to_end = as.numeric(day_to_end),
                yr_from_dm = as.numeric(yr_from_dm),
                yr_to_end = as.numeric(yr_to_end)) %>%
  bind_rows(demo_from_obs) %>%
  bind_rows(demo_from_pat)

other_fact_stack %>% filter(PATIENT_NUM == 70) %>% View

save(other_fact_stack, file="./data/DKD_heron_other_facts_prep.Rdata")


#RX at SCDF/SCBF(level4); DX at ICD concept(level5)
fact_stack <- other_fact_stack %>%
  dplyr::select(PATIENT_NUM,
                VARIABLE_CATEG,
                CONCEPT_CD,
                NVAL_NUM,
                START_DATE) %>%
  bind_rows(med_stack2 %>%
              dplyr::filter(CONCEPT_LEVEL==4) %>%
              mutate(CONCEPT_CD=paste0(CONCEPT_CD,"@",RX_MODIFIER)) %>%
              dplyr::select(PATIENT_NUM,
                            VARIABLE_CATEG,
                            CONCEPT_CD,
                            NVAL_NUM,
                            START_DATE)) %>%
  bind_rows(dx_stack2 %>%
              dplyr::filter(CONCEPT_LEVEL==5) %>%
              mutate(CONCEPT_CD=paste0(CONCEPT_CD,"@",DX_MODIFIER)) %>%
              dplyr::select(PATIENT_NUM,
                            VARIABLE_CATEG,
                            CONCEPT_CD,
                            NVAL_NUM,
                            START_DATE))

load("./data/pat_tbl.Rdata")

fact_stack %<>%
  left_join(pat_tbl %>% 
              dplyr::select(PATIENT_NUM,DM_DATE,END_DATE) %>% unique,
            by="PATIENT_NUM") %>%
  dplyr::mutate(day_from_dm=as.numeric(START_DATE-DM_DATE),
                day_to_end=as.numeric(END_DATE-START_DATE)) %>%
  dplyr::mutate(yr_from_dm=floor(day_from_dm/365.25),
                yr_to_end=floor(day_to_end/365.25))

View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
       dplyr::summarize(fact_cnt=n(),
                        pat_cnt=length(unique(PATIENT_NUM)),
                        cd_cnt=length(unique(CONCEPT_CD))))


## Preprocess 7 -- remove features with less than 5 nonDKD or DKD cases
overall_summary<-fact_stack %>%
  left_join((pat_tbl %>% dplyr::select(PATIENT_NUM,DKD_IND)),
            by="PATIENT_NUM") %>%
  group_by(VARIABLE_CATEG,CONCEPT_CD) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM)),
                   dkd_cnt = length(unique(PATIENT_NUM * DKD_IND))-1) %>%
  mutate(nondkd_cnt = pat_cnt,
         dkd_rt = dkd_cnt/pat_cnt,
         nondkd_rt = 1-dkd_cnt/pat_cnt) %>%
  mutate(odd_ratio_emp=dkd_rt/nondkd_rt) %>%
  arrange(desc(odd_ratio_emp))
  
cnt_before<-nrow(fact_stack)
fact_stack %<>% 
  semi_join((overall_summary %>% 
               dplyr::filter(dkd_cnt>=5 & nondkd_cnt>=5)),
            by="CONCEPT_CD")
cnt_after<-nrow(fact_stack)
#report
cat("Filtered out",cnt_before - cnt_after," facts due to normality assumption.\n")
cat("Filtered out",cd_before - cd_after," variables due to normality assumption.\n")
gc()

View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
       dplyr::summarize(fact_cnt=n(),
                        pat_cnt=length(unique(PATIENT_NUM)),
                        cd_cnt=length(unique(CONCEPT_CD))))

## Preprocess 8 -- frequency filter
pass_ratio<-0.01 #user-defined
pass_filter<-fact_stack %>%
  group_by(VARIABLE_CATEG) %>%
  do(mutate(.,fact_cnt_categ=nrow(.),
              pat_cnt_categ=length(unique(PATIENT_NUM)),
              cd_cnt_categ=length(unique(CONCEPT_CD)))) %>%
  group_by(VARIABLE_CATEG,CONCEPT_CD,
           pat_cnt_categ,fact_cnt_categ,cd_cnt_categ) %>%
  dplyr::summarize(pat_cnt = length(unique(PATIENT_NUM))) %>%
  dplyr::filter(pat_cnt >= min(pass_ratio*pat_cnt_categ,
                               fact_cnt_categ/cd_cnt_categ))

cnt_before<-nrow(fact_stack)
cd_before<-length(unique(fact_stack$CONCEPT_CD))
fact_stack %<>%
  semi_join(pass_filter, by=c("CONCEPT_CD"))
cnt_after<-nrow(fact_stack)
cd_after<-length(unique(fact_stack$CONCEPT_CD))
#report
cat("Filtered out",cnt_before - cnt_after," facts due to high sparsity.\n")
cat("Filtered out",cd_before - cd_after," variables due to high sparsity.\n")
gc()
# Filtered out 2,260,814 facts due to low frequency.
# Filtered out 18,017  variables due to low frequency.

# same-day averaging
fact_stack %<>%
  group_by(PATIENT_NUM,VARIABLE_CATEG,CONCEPT_CD,START_DATE,
           DM_DATE,END_DATE,day_from_dm,day_to_end,yr_from_dm) %>%
  dplyr::summarize(NVAL_NUM = mean(NVAL_NUM)) %>% ungroup

# View(fact_stack %>% group_by(VARIABLE_CATEG) %>%
#        dplyr::summarize(fact_cnt=n(),
#                         pat_cnt=length(unique(PATIENT_NUM)),
#                         cd_cnt=length(unique(CONCEPT_CD))))

save(fact_stack, file="./data/DKD_heron_facts_prep.Rdata")
