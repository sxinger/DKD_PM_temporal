##---------------------------helper functions--------------------------------------##

#########################
#all about initial setup#
#########################
#install (if needed) and require packages
require_libraries<-function(package_list){
  #install missing packages
  new_packages<-package_list[!(package_list %in% installed.packages()[,"Package"])]
  if(length(new_packages)>0){
    install.packages(new_packages)
  }
  
  for (lib in package_list) {
    library(lib, character.only=TRUE)
    cat("\n", lib, " loaded.", sep="")
  }
}


dgCS_write_to_svm<-function(sparse_mx,file){
  if(!("dgCMatrix" %in% class(sparse_mx))){
    stop("sparse_mx must be a dgCMatrix")
  }
  
  dt<-data.table(i=sparse_mx@i+1,
                 j=rep(1:sparse_mx@Dim[2],each=sparse_mx@Dim[1]),
                 x=sparse_mx@x)
  dt<-dt[(x!=0),.(i, jx = paste(j, x, sep = ":"))]
  dt<-dt[order(i),.(res = paste(jx, collapse = " ")), by = i][["res"]]
  
  out<-paste(seq_len(sparse_mx@Dim[1]),dt)
  writeLines(out, con = file)
  
  out <-h2o.uploadFile(file,parse_type="SVMLight")
  return(out)
}



#knit hooks configurations
#borrowed from Steve Simon (acknowledge)
# require("knitr")
# track_time<<-TRUE
# if(track_time){
#   knit_hooks$set(timer=function(before, options, envir){
#     if(before){
#       current_time <<- Sys.time()
#       m <- paste("Chunk ", options$label, " started at ", as.character(current_time), ".\n", sep="")
#       cat(m)
#       return(m)
#     }else{
#       elapsed_time <- Sys.time() - current_time
#       elapsed_message <- paste(round(elapsed_time, 1), attr(elapsed_time, "units"))
#       m <- paste("Chunk", options$label, "used", elapsed_message)
#       message(m)
#       
#       if(exists("timing_log")) {
#         timing_log <<- c(timing_log, m)
#       }else{
#         timing_log <<- m
#       }
#       
#       write.table(timing_log, file="timing_log.txt", row.names=FALSE, col.names=FALSE) 
#       return(m)
#     }})
# }


#set up local instance using sparklyr
#deid: an csv file with 
setup_spark_local<-function(){
  #install spark locally
  spark_install(version="2.2.0") #install to default directory
  if (nchar(Sys.getenv("SPARK_HOME")) < 1) {
    Sys.setenv(SPARK_HOME = file.path(spark_install_dir(),spark_installed_versions()$dir)) #point to the directory where spark was installed
  }
  
  #setup connection with oracle db
  oracle.jar <-file.path("/usr/lib/oracle/11.2/client64/lib","ojdbc6.jar")
  config <- spark_config()
  config$`sparklyr.shell.driver-class-path`<- oracle.jar
  
  #configure spark
  config$`sparklyr.shell.driver-memory` <- "16G"
  config$spark.memory.fraction <- 0.8 #utilize more memory per executor
  
  #establish local connection
  sc <- spark_connect(master = "local", config = config)
  out<-list(config=config, sc=sc)
  return(out)
}

#borrow from Steve's 'archive' function, with some customization
archive<- function(dt, note="No notes", idx="",write.in=F) {
  # get the name of the data table
  dt_name <- paste(deparse(substitute(dt)), idx, sep="")
  if (!exists("arc")) {
    arc <<-list(documentation=NULL)
  }
  arc[[dt_name]] <<- dt
  name_with_dimensions <- paste(c(dt_name, dim(dt)), collapse=" ")
  arc$documentation <<- rbind(arc$documentation, data.table(name=name_with_dimensions, note=note))
  items_in_archive <- length(names(arc))
  size_of_archive <- format(object.size(arc), units="auto")
  cat(dt_name, " added. Archive has",
      items_in_archive-1, "items and uses",
      size_of_archive, ".\n")
  print(arc$documentation)
  return(invisible(names(arc))) # invisible prevents the names from printing again.
  
  if(write.in){
    save(dt,file=paste0('./data/',dt_name,'.Rdata'))
  }
}

##########################
#simple utility functions#
##########################

#find mode;
get_mode<-function(x){
  ux<-unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#transpose data.table
t_dt<-function(dt){
  dt<-as.data.table(dt)[,id_row:=1:.N]
  dt<-dcast(melt(dt,id.vars = "id_row"), variable~id_row)
  dt[,id_row:=NULL]
  return(dt)
}

#find all high-correlated pair
all_highcorr_pair<-function(corr_dt,cutoff=0.9){
  results<-data.table(v1=character(0),v2=character(0),
                      corr=numeric(0),stringsAsFactors=F)
  diag(corr_dt)<-0
  maxval<-1
  while(maxval>=cutoff){
    maxval<-max(abs(corr_dt))
    max_idx<-which(corr_dt==maxval,arr.ind=T)[1,]
    results<-rbind(results, data.table(v1=rownames(corr_dt)[max_idx[1]],
                                       v2=colnames(corr_dt)[max_idx[2]],
                                       corr=maxval))
    corr_dt[max_idx[1],max_idx[2]]<-0
    corr_dt[max_idx[2],max_idx[1]]<-0
  }
  return(results)
}

#in-place replacement of old values by single new value
update_nto1<-function(dt,olds,new){
  for(j in seq_along(dt)){
    set(dt,i=which(dt[[j]] %in% olds),j=j,value=new)
  }
  return(dt)
}

#convert long skinny data.frame to wide sparse matrix
long_to_sparse_matrix<-function(df,id,variable,val,binary=FALSE){
  if(binary){
    x_sparse<-with(df,
                   Matrix::sparseMatrix(i=as.numeric(as.factor(get(id))),
                                j=as.numeric(as.factor(get(variable))),
                                x=1,
                                dimnames=list(levels(as.factor(get(id))),
                                              levels(as.factor(get(variable))))))
  }else{
    x_sparse<-with(df,
                   Matrix::sparseMatrix(i=as.numeric(as.factor(get(id))),
                                j=as.numeric(as.factor(get(variable))),
                                x=ifelse(is.na(get(val)),1,as.numeric(get(val))),
                                dimnames=list(levels(as.factor(get(id))),
                                              levels(as.factor(get(variable))))))
  }
  
  return(x_sparse)
}


#get coefficents and their ranks
get_coef_rank<-function(fit,s="lambda.1se",...){
  get_coef<-coef(fit,s=s)
  coef_dt<-data.table(Feature=row.names(get_coef), #same name as in xgb.var.importance
                      coef_val=get_coef[,1],
                      Gain=abs(get_coef[,1])) #same name as in xgb.var.importance
  coef_dt[,coef_rank:=rank(-Gain)]
  coef_dt<-coef_dt[(coef_val!=0)]; gc()
  return(coef_dt[order(coef_rank)])
}


#####################
#data transformation#
#model tuning       #
#model validation   #
#####################

# multiple ggplots_lst sharing legends
require_libraries(c( "grid"
                    ,"gridExtra"))

grid_arrange_shared_legend<-function(plots_lst = list(),
                                     ncol = ncol,
                                     nrow = nrow,
                                     position = c("bottom", "right"),
                                     plot=T) {
  if(!is.list(plots_lst)){
    stop("Plots should be supplied in a list!")
  }
  position<-match.arg(position)
  g<-ggplotGrob(plots_lst[[1]] + 
                  theme(legend.position = position))$grobs
  legend<-g[[which(sapply(g,function(x) x$name) == "guide-box")]]
  lheight<-sum(legend$height)
  lwidth<-sum(legend$width)
  
  # nullify legends for each individual plot
  gl<-lapply(plots_lst,
             function(x) x + theme(legend.position = "none"))
  gl<-c(gl,
        ncol = ncol,
        nrow = nrow)
    
  combined<-switch(
    position,
    "bottom" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight)
    ),
    "right" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 2,
      widths = unit.c(unit(1, "npc") - lwidth, lwidth)
    )
  )
  
  if(plot){
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
  }else{
    return(combined)
  }

}

#lay out mutiple plots on one canvas
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#plot distribution over time for Sepsis project
plot_dist_over_time<-function(conn,event_tbl,start_event="",start_ref="",
                              critical_event,time_col,x_trunc,x_breaks){
  sql_query<-paste0("with event_order as (
                     select ENCOUNTER_NUM, EVENT_TYPE, EVENT_SUBTYPE, EVENT_SUBTYPE2, ",time_col,", 
                            rank() over (partition by ENCOUNTER_NUM order by ",time_col," asc) rn
                     from ",event_tbl,")
                         ,start_event as (
                     select ENCOUNTER_NUM, rn 
                     from event_order
                     where ",start_ref,")
                     
                     select x.EVENT_TYPE, x.EVENT_SUBTYPE, x.EVENT_SUBTYPE2, x.",time_col,",
                            case when x.EVENT_TYPE like 'Sepsis%' then 'Alert' else 'Treatment' end as Alert_or_Treat
                     from event_order x
                     where exists (select 1 from start_event y
                                   where y.ENCOUNTER_NUM = x.ENCOUNTER_NUM and
                                         y.rn <= x.rn) and
                           x.EVENT_SUBTYPE not like '6-hour%'")
  
  dat<-dbGetQuery(conn,sql_query) %>%
    dplyr::filter(get(time_col) <= x_trunc) %>%
    dplyr::mutate(time_1digit = round(get(time_col),1))
  
  ggplot(dat,aes(time_1digit))+
    geom_density(aes(fill=factor(get(critical_event)),
                     linetype=factor(ALERT_OR_TREAT)),alpha=0.5)+
    geom_rug(sides="b")+
    labs(title="Density plot", 
         subtitle=paste0("Hours since triage after ",start_event," is observed"),
         x="Hours since Triage",
         fill="Critical Events",
         linetype="Sepsis Alert or Treatment Steps")+
    scale_x_continuous(breaks=x_breaks)
}

#plot survival ROC
plot_survivalROC <- function(pred_df,t) {
  res <- with(pred_df,
              survivalROC(Stime        = time,
                          status       = status,
                          marker       = lp,
                          predict.time = t,
                          method       = "KM"))       # KM method without smoothing
  # res_smooth<-with(pred_df,
  #                  risksetROC(Stime        = time,
  #                             status       = status,
  #                             marker       = lp,
  #                             predict.time = t,
  #                             plot=F))                    # KM method with smoothing
  
  ## Plot ROCs
  plot(TP ~ FP, res, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC))
  # lines(TP ~ FP, res_smooth, type="l",lty=3)
  abline(a = 0, b = 1, lty = 2)
}

##grid search for alpha in Elastic Net
glmnet_alpha_cv<-function(x,y,alpha_seq,family,type.measure="deviance",s="lambda.min",pred_type="link",
                          maximize=F, nfolds_alpha=5,nfolds_lambda=5,parallel=T,cores,verb=T)
{
  y<-drop(y)
  N<-nrow(x)
  fold_alpha<-sample(seq_len(nfolds_alpha),N,replace=T)
  
  pred<-matrix(NA,nrow=N,ncol=length(alpha_seq))
  cverr<-c()
  lam_alpha<-matrix(NA,nrow=nfolds_alpha,ncol=length(alpha_seq),
                    dimnames=list(paste0("fold_",seq_len(nfolds_alpha)),
                                  paste0("alpha_",alpha_seq)))
  lam_alpha_cvm<-lam_alpha
  lam_alpha_k<-lam_alpha_cvm
  
  require_libraries(c("parallel","doParallel"))
  cl <- makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  for(j in seq_along(alpha_seq)){
    start_j<-Sys.time()
    if(verb){
      cat('alpha =',alpha_seq[j],'\n')
    }
    for(i in 1:nfolds_alpha){
      start_i<-Sys.time()
      if(verb){
        cat('...alpha fold:',i,'\n')
      }
      which_fold = (fold_alpha==i)
      if(family=="cox"){
        cv_glmnet<-cv.glmnet(x[!which_fold,],y[!which_fold,],family=family,type.measure=type.measure,
                             parallel=parallel,nfolds = nfolds_lambda,alpha=alpha_seq[j])
      }else{
        cv_glmnet<-cv.glmnet(x[!which_fold,],y[!which_fold],family=family,type.measure=type.measure,
                             parallel=parallel,nfolds = nfolds_lambda,alpha=alpha_seq[j])
      }
      
      if(s=='lambda.min'){
        lam_alpha[i,j]<-cv_glmnet$lambda.min
        lam_alpha_cvm[i,j]<-cv_glmnet$cvm[which.min(cv_glmnet$lambda)]
      }else{
        lam_alpha[i,j]<-cv_glmnet$lambda.1se
        lam_alpha_cvm[i,j]<-cv_glmnet$cvm[cv_glmnet$lambda==cv_glmnet$lambda.1se]
      }
      lam_alpha_k[i,j]<-nrow(predict(cv_glmnet,s=s,type="nonzero"))-1
      pred[which_fold,j]<-predict(cv_glmnet,
                                  newx=x[which_fold,],
                                  s=lam_alpha[i,j],
                                  type=pred_type)
      if(verb){
        lapse_i<-Sys.time()-start_i
        cat('...finish tune lambda for alpha fold:',i,'in',lapse_i,units(lapse_i),'\n')
      }
    }
    if(family=="cox"){
      calibr<-coxph(Surv(time,status)~offset(lp),
                    data=data.frame(time=y[,1],status=y[,2],lp=pred[,j]))
      cverr<-calibr$loglik[1]
    }else{
      calibr<-glm(y~offset(lp),
                  data=data.frame(y=y,lp=pred[,j]),
                  family=family)
      cverr<-2*(median(lam_alpha_k[,j]) - calibr$null.deviance)
    }
    
    if(verb){
      lapse_j<-Sys.time()-start_j
      cat('end alpha =',alpha_seq[j],'in',lapse_j,units(lapse_j),'\n')
    }
    gc()
  }
  
  if(maximize){
    alpha_opt<-alpha_seq[which.max(cverr)]
    lambda_opt<-lam_alpha[which.max(lam_alpha_cvm[,which.max(cverr)]),which.max(cverr)]
  }else{
    alpha_opt<-alpha_seq[which.min(cverr)]
    lambda_opt<-lam_alpha[which.min(lam_alpha_cvm[,which.min(cverr)]),which.min(cverr)]
  } 
  
  outlist<-data.frame(alpha_opt=alpha_opt,
                      lambda_opt=lambda_opt)
  return(outlist)
  
  stopCluster(cl)
}

xgb_tune_hyper_param<-function(dtrain,params_grid,eval_metric,objective,nfold,maximize=T,booster="gbtree",
                               nrounds=2000,early_stopping_rounds=200,print_every_n=10,verb=T){
  
  bst_grid<-c()
  metric_name<-paste0("test_", eval_metric,"_mean")
  metric_sd_name<-paste0("test_", eval_metric,"_std")
  
  for(i in seq_len(dim(params_grid)[1])){
    start_i<-Sys.time()
    param<-as.list(params_grid[i,])
    
    bst <- xgb.cv(param,dtrain,nrounds=nrounds,nfold=nfold,objective = objective, metrics = eval_metric,booster=booster,
                  maximize=maximize,early_stopping_rounds= early_stopping_rounds, print_every_n = print_every_n)
    if(maximize){
      bst_grid<-rbind(bst_grid, cbind(params_grid[i,],
                                      bst$evaluation_log[which.max(bst$evaluation_log[[metric_name]]),]))
    }else{
      bst_grid<-rbind(bst_grid, cbind(params_grid[i,],
                                      bst$evaluation_log[which.min(bst$evaluation_log[[metric_name]]),]))
    }
    
    if(verb){
      cat("finished train case:",paste(colnames(params_grid),param,collapse="; "),
          "based on",eval_metric,"in",Sys.time()-start_i,units(Sys.time()-start_i),".\n")
      start_i<-Sys.time()
    }
  }
  return(bst_grid)
}

##====form the data into survival-object format====##
#fixed intervals:every k days, [k,2k)
#input data set requirements:
# - id: identification key (should be the same between start_end and all_obs)
# - start: start point
# - end: case/censor point
# - date: stack all observation dates
# - status: binary, 0-censor(could be status for cases prior to their case point); 1-case
# - covar: covariates(predictors) that change over time

as_PWSurvObj<-function(start_end,all_obs,id="",start="",end="",date="",status="",covar){
  start_end<-as.data.table(start_end); gc() #convert to data.table
  all_obs<-as.data.table(all_obs)[,date_type:=date]; gc() #convert to data.table
  setnames(start_end,status,"status") 
  setnames(all_obs,date,"date")
  
  start_end_long<-melt(start_end,id.vars=id,measure.vars=c(start,end),
                       variable.name="date_type",value.name='date')
  date_dt<-rbind(start_end_long,all_obs[,c(id,"date_type","date",covar),with=F],fill=T)
  date_dt<-merge(date_dt,start_end[,c(id,end,"status"),with=F],
                 by.x=c(id,"date"),by.y=c(id,end),all=T)
  date_dt[is.na(status),status:=0]
  date_dt<-update_nto1(date_dt,NaN,NA)
  gc()
  
  return(date_dt[order(get(id),date)])  
}

#input data set requirements:
# - id: identification key
# - date: stack all observation dates, start and end dates
# - status: binary, 0-censor(could be status for cases prior to their case point); 1-case
# - interval_size: size of observation intervals (fixed)
# - covar: covariates(predictors) that change over time
# - fun_lst: a list of function for aggregating covar
# - fun_str: the corresponding names of functions in fun_lst (strings)
create_fixed_interval<-function(dat,id="",date="date",status="status",interval_size,covar,
                                fun_lst,fun_str,impute_type=c("forward","backward","zero")){
  dat_cp<-data.table(dat); gc() #convert to data.table
  setnames(dat_cp,c(id,date,status),c("id","date","status"))
  
  dat_time<-dat_cp[,c("id","date","status"),with=F]
  setkeyv(dat_time,c("id","date","status"))
  dat_time<-unique(dat_time)
  gc()
  
  dat_time<-dat_time[order(id,date)]
  dat_time[,ffill:=date[1],by=id]
  dat_time[,time1:=as.numeric(date-ffill)] # in days
  dat_time[,c("time2","time3"):=list(floor(time1/interval_size)*interval_size,
                                     ceiling(time1/interval_size)*interval_size)]
  dat_cp<-merge(dat_cp,dat_time,by=c("id","date","status"))
  
  dat_pre<-dat_cp[,.(time1=time1[.N],time3=time3[.N],status=status[.N]),by=.(id,time2)]
  dat_pre<-dat_pre[,c("id","status","time2","time3","time1"),with=F]
  
  dat_other<-dat_cp[,mapply(do.call,fun_lst,lapply(.SD,list)),
                 .SDcols=covar,by=.(id,time2)]
  other_var_agg<-paste(paste0(fun_str,"_"),covar,sep="")
  dat_other[,variable:=other_var_agg]
  dat_other<-dcast(dat_other,...~variable,value.var="V1")
  dat1<-merge(dat_pre,dat_other,by=c("id","time2"))
  
  dat1_reconst<-dat1[,.(nbr_interval=time3[.N]/interval_size),by=id]
  dat1_reconst<-dat1_reconst[(rep(1:.N,nbr_interval))]
  dat1_reconst[,time3:=(1:.N)*interval_size,by=id]
  dat1_reconst[,time2:=shift(time3,n=1L,fill=0),by=id]
  dat1_reconst[,time2:=as.numeric(time2)]
  
  dat2<-merge(dat1_reconst,dat1,by=c("id","time2","time3"),all.x=T)
  vari<-colnames(dat2[,-c("id","time2","time3","nbr_interval"),with=F])
  
  if(impute_type=="zero"){
    dat2[,(vari):=lapply(.SD,function(x) ifelse(is.na(x),0,x)),.SDcols=vari]
  }else if(impute_type=="forward"){
    dat2[,(vari):=lapply(.SD,function(x) na.locf(x,na.rm=F)),.SDcols=vari,by=.(id)]
  }else if(impute_type=="backward"){
    dat2[,(vari):=lapply(.SD, function(x) na.locf(x,na.rm=F,fromLast=T)),.SDcols=vari,by=.(id)]
  }
  gc()
  
  setnames(dat2,c("id"),c(id))
  return(dat2)
}

#input data set requirements:
# - id: identification key
# - date: stack all observation dates, start and end dates
# - status: binary, 0-censor(could be status for cases prior to their case point); 1-case
# - interval_size: size of observation intervals (fixed)
# - covar: covariates(predictors) that change over time
# - fun_lst: a list of function for aggregating covar
# - fun_str: the corresponding names of functions in fun_lst (strings)
create_flex_interval<-function(dat,id="",date="date",status="status",covar,
                               impute_type=c("forward","backward","zero")){
  dat_cp<-data.table(as.data.table(dat)[,c(id,date,status,covar),with=F]); gc() #keep neccessary columns and make a copy
  setnames(dat_cp,c(id,date,status),c("id","date","status")) #rename
  
  #remove replicates(replicates-like, exactly same observations(or NA) occur on the same day)
  setkeyv(dat_cp,colnames(dat_cp))
  dat_cp<-unique(dat_cp); gc()
  
  dat_cp<-dat_cp[order(id,date)]; gc()
  dat_cp[,ffill:=date[1],by=id]
  
  dat_cp[,time2:=as.numeric(date-ffill)] # in days
  dat_cp[,time1:=shift(time2,n=1L,type="lag",fill=0),by=.(id)] # get end point of each observation interval
  
  if(impute_type=="zero"){
    dat_cp[,(covar):=lapply(.SD,function(x) ifelse(is.na(x),0,x)),.SDcols=covar]
  }else if(impute_type=="forward"){
    dat_cp[,(covar):=lapply(.SD,function(x) na.locf(x,na.rm=F)),.SDcols=covar,by=.(id)]
    dat_cp[,(covar):=lapply(.SD,function(x) na.locf(x,fromLast=T)),.SDcols=covar,by=.(id)]
  }else if(impute_type=="backward"){
    dat_cp[,(covar):=lapply(.SD, function(x) na.locf(x,na.rm=F,fromLast=T)),.SDcols=covar]
    dat_cp[,(covar):=lapply(.SD,function(x) na.locf(x,fromLast=F)),.SDcols=covar,by=.(id)]
  }
  gc()
  
  setnames(dat_cp,c("id","date","status"),c(id,date,status))
  dat_cp<-dat_cp[(time1+time2>0),c(id,"time1","time2",status,covar),with=F]; gc() #reorder columns for better look
  return(dat_cp)
}

#partial effect table for 'black box' models
partial_effect<-function(fit,x,nrs=30,rsp=0.5,pred_grid=NULL,pred_var,num_interval=20){
  if(is.null(pred_grid)){
    pred_grid<-partial(object=fit,
                       grid.resolution = min(nlevels(as.factor(x[,pred_var])),num_interval),
                       pred.var=pred_var,
                       train=x,
                       prob = T,
                       progress="text") %>%
      dplyr::filter(get(pred_var)!=0) %>%
      dplyr::select(pred_var)
  }else{
    pred_grid<-data.frame(pred=pred_grid)
    colnames(pred_grid)[1]<-pred_var
  }
  
  get_peff<-c()
  for(i in 1:nrs){
    cat('resample:',i,'for evaluating partial effect for',pred_var,'\n')
    resamplex<-x[(base::sample(c(TRUE,FALSE),
                               nrow(x),
                               c(rsp,1-rsp),replace=T)),]
    
    part_effect<-partial(object=fit,
                         pred.var=pred_var,
                         pred.grid = pred_grid,
                         train=resamplex,
                         prob = T,
                         progress="text")
    get_peff<-get_peff %>% 
      bind_rows(part_effect)
  }

  return(get_peff)
}

