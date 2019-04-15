rm(list=ls())

source("./R/util.R")

require_libraries(c("DBI",
                    "ROracle",
                    "dplyr",
                    "tidyr",
                    "magrittr"))

#send patient_num, encounter_num over to nheron and retrieve current patient_num, encounter_num
pat_num<-readRDS("./data2/pat_episode2.rda") %>%
  dplyr::select(PATIENT_NUM) %>% unique

config<-read.csv("./config_bh.csv")
conn<-connect_to_oci(config)

dbWriteTable(conn,"DKD_PAT_NUM_CALAMUS",pat_num,overwrite=T)


##..............work on sql oracle...............##


## Q1-how many T1DM and CFRD? necessary to re-run the whole experiment?
t1dm<-dbGetQuery(conn,"select * from T1DM_PAT_CALAMUS")
saveRDS(t1dm,file="./data2/pat_T1DM.rda")


