########utility functions##########
# handle temporal data
# require:
# - fact_stack
# - pat_episode
# - ep_unit: episode units in days


get_most_recent<-function(pat_episode,fact_stack,ep_unit){
  #### non-temporal model -- most recent value
  x_add<-pat_episode %>%
    #each patient has multiple episodes
    dplyr::select(PATIENT_NUM,episode) %>% 
    dplyr::group_by(PATIENT_NUM) %>% 
    #take their last available episode
    top_n(n=1,wt=episode) %>% 
    #use onset episode indicator as additional feature
    dplyr::mutate(CONCEPT_CD=paste0("ep_",episode),
                  NVAL_NUM=1) %>% 
    ungroup %>% unique
  
  X_long<-fact_stack %>% 
    inner_join(x_add %>% dplyr::select(PATIENT_NUM,episode) %>% unique,
               by="PATIENT_NUM") %>%
    #collect all historical values at least 1-episode prior to target
    filter(day_from_dm < episode*ep_unit) %>% 
    dplyr::group_by(PATIENT_NUM,CONCEPT_CD) %>% 
    #use most recent value
    dplyr::top_n(n=-1,wt=day_to_end) %>% ungroup %>% 
    dplyr::select(PATIENT_NUM, CONCEPT_CD, NVAL_NUM) %>% 
    #add years since DM onset
    bind_rows(x_add %>% dplyr::select(-episode)) %>% 
    #sort by patient_num
    dplyr::arrange(PATIENT_NUM) %>% 
    unique #remove duplicates
  
  #for alignment with target
  x_pn<-x_add %>% dplyr::select(PATIENT_NUM) %>%
    unique %>% arrange(PATIENT_NUM)
  
  #convert to sparse matrix
  X_wide <- X_long %>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM",
                          variable="CONCEPT_CD",
                          val="NVAL_NUM")
  rm(X_long); gc()
  
  #collect target
  y_pn<-pat_episode %>%
    semi_join(x_pn,by="PATIENT_NUM") %>% #
    dplyr::select(PATIENT_NUM,DKD_IND_additive) %>%
    unique %>% group_by(PATIENT_NUM) %>% 
    top_n(n=1,wt=DKD_IND_additive) %>% ungroup
  
  #alignment check
  align_row<-all((row.names(X_wide)==y_pn$PATIENT_NUM)) # yes
  
  if(!align_row) {
    stop("rows for convariate matrix and target don't align!")
  }
  
  Xy_all<-list(X_wide=X_wide,
               y_pn=y_pn)
  return(Xy_all)
}


get_stack_temporal<-function(pat_episode,X_long){
  #pivot to sparse matrix
  X_long %<>%
    semi_join(pat_episode, by="PATIENT_NUM") %>%
    unite("VARIABLE_ep",c("VARIABLE","episode_x")) %>%
    dplyr::select(PATIENT_NUM, VARIABLE_ep, NVAL_NUM) %>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM",
                          variable="VARIABLE_ep",
                          val="NVAL_NUM")
  
  X_idx<-data.frame(PATIENT_NUM = as.integer(row.names(X_long)),
                    stringsAsFactors = F)
  
  #collect target
  y_long<-pat_episode %>%
    semi_join(X_idx,by="PATIENT_NUM") %>% #
    dplyr::select(PATIENT_NUM,episode,DKD_IND_additive) %>%
    group_by(PATIENT_NUM) %>% 
    top_n(n=1,wt=episode) %>% ungroup %>% 
    dplyr::select(PATIENT_NUM,DKD_IND_additive) %>%
    group_by(PATIENT_NUM) %>% 
    top_n(n=1,wt=DKD_IND_additive) %>% ungroup %>%
    unique
  
  #alignment check
  align_row<-all((row.names(X_long)==y_long$PATIENT_NUM)) # yes
  
  if(!align_row) {
    stop("rows for convariate matrix and target don't align!")
  }
  
  Xy_all<-list(X_ep = X_long,
               y_ep = y_long)
  
  return(Xy_all)
}


get_dsurv_temporal<-function(pat_episode,X_long){
  #pivot to sparse matrix
  X_long %<>%
    left_join(pat_episode %>% group_by(PATIENT_NUM, episode) %>% 
                top_n(n=1,wt=DKD_IND_additive) %>% ungroup %>%
                dplyr::select(PATIENT_NUM,episode) %>% unique,
              by = "PATIENT_NUM") %>%
    dplyr::filter(episode_x < episode) %>% 
    unite("VARIABLE_ep",c("VARIABLE","episode_x")) %>%
    arrange(PATIENT_NUM, episode) %>%
    unite("PATIENT_NUM_ep",c("PATIENT_NUM","episode")) %>%
    dplyr::select(PATIENT_NUM_ep, VARIABLE_ep, NVAL_NUM) %>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM_ep",
                          variable="VARIABLE_ep",
                          val="NVAL_NUM")
  
  X_idx<-data.frame(PATIENT_NUM_ep = row.names(X_long),
                    stringsAsFactors = F)
  
  #collect target
  y_long<-pat_episode %>%
    dplyr::select(PATIENT_NUM,episode, DKD_IND_additive) %>% unique %>%
    group_by(PATIENT_NUM, episode) %>% 
    top_n(n=1,wt=DKD_IND_additive) %>% ungroup %>%
    unite("PATIENT_NUM_ep",c("PATIENT_NUM","episode")) %>%
    semi_join(X_idx,by="PATIENT_NUM_ep") %>%
    arrange(PATIENT_NUM_ep)
  
  #alignment check
  align_row<-all((row.names(X_long)==y_long$PATIENT_NUM_ep)) # yes
  
  if(!align_row) {
    stop("rows for convariate matrix and target don't align!")
  }
  
  Xy_all<-list(X_ep = X_long,
               y_ep = y_long)
  
  return(Xy_all)
}



library(Matrix)
get_grptemporal<-function(pat_episode,X_long,verb=TRUE){
  # transformation should only be derived from training
  X_long_tr<-X_long %>% 
    semi_join(pat_episode %>% filter(!(part_idx %in% c("H","5"))),
              by=c("PATIENT_NUM"))
  
  # feature-wise standardization (output)
  center_scale<-X_long_tr %>%
    group_by(VARIABLE) %>%
    dplyr::summarize(NVAL_mean=mean(NVAL_NUM,na.rm=T),
                     NVAL_sd=ifelse(length(NVAL_NUM)==1,0,sd(NVAL_NUM,na.rm=T))) %>%
    filter(NVAL_sd > 1e-5) #filter out really low variance
  
  # standardization 
  X_long_tr %<>%
    inner_join(center_scale,by=c("VARIABLE")) %>%
    dplyr::mutate(NVAL_std = (NVAL_NUM - NVAL_mean)/NVAL_sd)
  
  # time-wise orthognalization
  n<-length(unique(X_long_tr$PATIENT_NUM))
  G<-unique(X_long_tr$VARIABLE)
  J<-length(G)
  Tr<-vector("list", J)
  Xstd_long_tr<-c()
  for(j in seq_along(numeric(J))){
    if(verb){
      start_j<-Sys.time()
      cat("transform variable",G[j],"across time windows.\n")
    }

    ind<-X_long_tr %>% filter(VARIABLE==G[j]) %>%
      dplyr::select(episode_x) %>% 
      arrange(episode_x) %>%
      unique %>% unlist
    
    if (length(ind)==0) next
    X<-X_long_tr %>% filter(VARIABLE==G[j]) %>%
      dplyr::select(PATIENT_NUM,VARIABLE,episode_x,NVAL_std)
    
    if(length(ind)==1){
      Xstd_long_tr %<>% 
        bind_rows(X %>% dplyr::select(-episode_x))
    }else{
      X %<>% 
        unite("VARIABLE_ep",c("VARIABLE","episode_x")) %>%
        long_to_sparse_matrix(.,
                              id="PATIENT_NUM",
                              variable="VARIABLE_ep",
                              val="NVAL_std")
      
      SVD<-svd(X, nu=0)
      r <- which(SVD$d > 1e-2)
      Tr[[j]] <- sweep(SVD$v[,r,drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
      XX<-X%*%Tr[[j]]
      colnames(XX)<-paste0(G[j],"_Comp",r)
      nz <- !apply(XX==0,2,all)
      XX <- XX[, nz, drop=FALSE]
      
      Xstd_long_tr %<>%
        bind_rows(data.frame(PATIENT_NUM = as.numeric(rep(row.names(XX),length(r))),
                             VARIABLE    = rep(colnames(XX),each=nrow(XX)),
                             NVAL_std    = XX@x,
                             stringsAsFactors = F))
    }

    if(verb){
      lapse_j<-Sys.time()-start_j
      cat("finish transform variable",G[j],"across time windows in",lapse_j,units(lapse_j),".\n")
    }
  }
  
  Xstd_long_tr %<>%
    long_to_sparse_matrix(.,
                          id="PATIENT_NUM",
                          variable="VARIABLE",
                          val="NVAL_std")
  names(Tr)<-G

  #for alignment with target
  X_idx<-data.frame(PATIENT_NUM = unique(X_long$PATIENT_NUM),
                    stringsAsFactors = F) %>%
    arrange(PATIENT_NUM)
  
  #collect target
  y_long<-pat_episode %>%
    semi_join(X_idx,by="PATIENT_NUM") %>% #
    dplyr::select(PATIENT_NUM,episode,DKD_IND_additive) %>%
    group_by(PATIENT_NUM) %>% 
    top_n(n=1,wt=episode) %>% ungroup %>% 
    dplyr::select(PATIENT_NUM,DKD_IND_additive) %>%
    group_by(PATIENT_NUM) %>% 
    top_n(n=1,wt=DKD_IND_additive) %>% ungroup %>%
    unique %>% arrange(PATIENT_NUM)
  
  #alignment check
  align_row<-all((X_idx$PATIENT_NUM==y_long$PATIENT_NUM)) # yes
  
  if(!align_row) {
    stop("rows for convariate matrix and target don't align!")
  }
  
  Xy_all<-list(X_ep = Xstd_long_tr,
               y_ep = y_long,
               Tr_std = center_scale,
               Tr_ort = Tr)
  
  return(Xy_all)
}
