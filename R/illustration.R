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
