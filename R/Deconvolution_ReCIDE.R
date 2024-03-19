
ReCIDE_deconvolution <- function(Sig_list,EXP_df,Method = 'DWLS',n_cores=5,n_celltype=0){
  
  EXP=list()
  for (i in 1:ncol(EXP_df)) {
    EXP[[i]]=as.numeric(EXP_df[,i])
    names(EXP[[i]])=row.names(EXP_df)
  }
  names(EXP)=colnames(EXP_df)
  
  
  
  
  if(length(Sig_list)<1.5*length(EXP)){
    
    fun_DWLS_in<-function(i){
      bulk<-EXP[[i]]
      sep_ref.list<-Sig_list
      
      sep_solDWLS_inside<-list()
      for (j in 1:length(sep_ref.list)) {
        ##dwls??????dataframe?????????matrix
        ref1<-as.matrix(sep_ref.list[[j]])
        
        # batch_output<-cosine_screen_HighToLow(ref1,bulk)
        batch_output<-cosine_screen_LowToHigh(ref1,bulk)
        # batch_output<-cosine_screen_GA(ref1,bulk)
        
        tr <- batch_output[[2]]
        # tr<-trimData(ref1,bulk)
        sep_solDWLS_inside[[j]]<-try(solveDampenedWLS(tr[[1]],tr[[2]]), TRUE)
        
        
        names(sep_solDWLS_inside)[j]<-names(sep_ref.list)[j]
        
        
      }
      return(sep_solDWLS_inside)
    }
    sep_solDWLS=pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores=n_cores)
    names(sep_solDWLS)=names(EXP)
    
  }else{
    sep_solDWLS=list()
    
    
    for (f in 1:length(EXP)) {
      bulk<-EXP[[f]]
      sep_ref.list<-Sig_list
      
      DWLS_output<-list()
      
      
      fun_DWLS_in2<-function(j){
        ref1<-as.matrix(sep_ref.list[[j]])
        
        batch_output<-cosine_screen_LowToHigh(ref1,bulk)
        
        tr <- batch_output[[2]]
        sep_solDWLS_j<-try(solveDampenedWLS(tr[[1]],tr[[2]]), TRUE)
        
        return(sep_solDWLS_j)
        
      }
      sep_solDWLS_j=mclapply(1:length(sep_ref.list),fun_DWLS_in2,mc.cores=n_cores)
      names(sep_solDWLS_j)=names(sep_ref.list)
      
      sep_solDWLS[[f]]=sep_solDWLS_j
    }
    names(sep_solDWLS)=names(EXP)
    
  }
  
  
  for (j in 1:length(sep_solDWLS)) {
    for(k in length(sep_solDWLS[[j]]):1){
      if(length(sep_solDWLS[[j]][[k]])<n_celltype){sep_solDWLS[[j]][[k]]<-NULL}
    }
  }
  
  # source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
  prd_df=ReCIDE_PCA(sep_solDWLS,method=Method)
  
  prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
  
  
  deconvolution_output=list(sep_solDWLS,prd_df)
  names(deconvolution_output)=c('results_list','results_final_df')
  return(deconvolution_output)
}





ReCIDE_PCA=function(results_deconvolution,method='DWLS'){
  
  for (i in length(results_deconvolution):1) {
    for (j in length(results_deconvolution[[i]]):1) {
      if(class(results_deconvolution[[i]][[j]])[1]=="try-error"){results_deconvolution[[i]][[j]]<-NULL}
    }
  }
  
  if(method=='FARDEEP'){
    for (i in 1:length(results_deconvolution)) {
      for (j in 1:length(results_deconvolution[[i]])) {
        results_deconvolution[[i]][[j]]<-results_deconvolution[[i]][[j]][["relative.beta"]][1,]
      }}
    prd<-results_deconvolution
  }
  ##DWLS????????????
  # # prd<-sep_solOLS
  if(method=='DWLS'){
    for (i in length(results_deconvolution):1) {
      for (j in length(results_deconvolution[[i]]):1) {
        if(class(results_deconvolution[[i]][[j]])=="try-error"){results_deconvolution[[i]][[j]]<-NULL}
      }
    }
    prd<-results_deconvolution
  }
  
  
  # cibersort????????????
  if(method=='CIBERSORT'){
    for (i in 1:length(results_deconvolution)) {
      for (j in 1:length(results_deconvolution[[i]])) {
        data_in=results_deconvolution[[i]][[j]]
        data_in<-as.data.frame(data_in[1,-((ncol(data_in)-2):ncol(data_in))])
        data_in=as.numeric(data_in)
        names(data_in)=colnames(results_deconvolution[[i]][[j]])[1:(ncol(results_deconvolution[[i]][[j]])-3)]
        results_deconvolution[[i]][[j]]=data_in
      }
    }
    prd<-results_deconvolution
  }
  
  #prd??????????????????df??????
  prd_df<-list()
  for(i in 1:length(prd)){
    
    prd[[i]][[1]]<-as.data.frame(prd[[i]][[1]])
    df_merge<-prd[[i]][[1]]
    colnames(df_merge)[1]<-names(prd[[i]])[1]
    
    for(j in 2:length(prd[[i]])){
      
      prd[[i]][[j]]<-as.data.frame(prd[[i]][[j]])
      df_merge<-merge(df_merge, prd[[i]][[j]], by = "row.names", all = TRUE)
      row.names(df_merge)<-df_merge[,1]
      df_merge<-df_merge[,-1]
      colnames(df_merge)[j]<-names(prd[[i]])[j]
    }
    
    prd_df[[i]]<-as.data.frame(lapply(df_merge,as.numeric))
    row.names(prd_df[[i]])<-row.names(df_merge)
    names(prd_df)[i]<-names(prd)[i]
  }
  
  
  
  #prd_df<-prd_dfFB
  prd_dfFB<-prd_df
  for (j in 1:length(prd_df)) {
    for(i in 1:length(prd_df[[j]][,1])){
      ###0.75?????????????????????
      if(table(is.na(prd_df[[j]][i,]))['FALSE']<0.75*ncol(prd_df[[j]])){
        prd_df[[j]][i,][is.na(prd_df[[j]][i,])]<-0}
    }
  }
  ####?????????PCA????????????????????????PCA??????
  ##PCA
  ####PCA?????? 
  
  #################
  # for (n in 1:length(prd_df)) 
  func_pca<- function(n){
    
    data_train <- t(prd_df[[n]])
    
    data_train[is.na(data_train)]=0
    
    for (i in ncol(data_train):1) {
      #??????????????????0???????????????????????????????????????
      if((nrow(as.data.frame(data_train[data_train[,i]<=0,]))>0.5*nrow(data_train))){
        data_train=data_train[,-i]}
    }
    data_train=as.data.frame(data_train)
    if(ncol(data_train)>=3){
      library(gmodels)
      #????????????????????????????????????????????????
      PCA.OUT=fast.prcomp(data_train, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL)
      
      imp=PCA.OUT[["sdev"]]/sum(PCA.OUT[["sdev"]])
      
      
      pc_sure=2
      
      pca_n=PCA.OUT$x[,1:pc_sure]
      
      pca_n[,1]=scale(pca_n[,1])
      pca_n[,2]=scale(pca_n[,2])
      
      ##pc???????????????10%
      ##???????????????70%
      ##???PC?????????10
      ##????????????1???mclust
      #####################################
      #####################################
      model <- Mclust(pca_n)
      
      
      # & (model_class[1,2]/nrow(pca_n)>0.5)
      ###????????????PCA??????????????????
      if(model$G !=1){
        model_class=as.data.frame(table(model$classification))
        model_class=model_class[order(model_class[,'Freq'],decreasing = TRUE),]
        
        if(model_class[1,2]>5){
          if(model_class[1,2]==model_class[2,2]){
            ##?????????????????????????????????
            
            model_class=model_class[model_class[,2]==model_class[1,2],]
            
            use_ind=model$classification[model$classification %in% model_class[,1]]
            data_train=data_train[names(use_ind),]
            
            data_train_list=split(data_train,use_ind) 
            
            kappa_value=c()
            for (ka in 1:length(data_train_list)) {
              data_in=as.data.frame(t(data_train_list[[ka]]))
              kappa_value[ka]=kappa(as.matrix(data_in))
            }
            
            # kappa<100????????????????????????????????????
            # 100<=kappa<=1000??????????????????????????????????????????????????????
            # ???kappa>1000????????????????????????????????????
            
            patient_names=row.names(data_train_list[[which.max(kappa_value)]])
          }else{
            ##?????????????????????????????????
            patient_names=names(model$classification[model$classification==model_class[1,1]])
          }
        }else{
          ##??????else??????????????????5
          patient_names=colnames(prd_df[[n]])
        }
        
      }else{
        ##?????????????????????1???????????????
        patient_names=colnames(prd_df[[n]])}
    }
    
    return_data=prd_df[[n]][,patient_names]
    return(return_data)
  }
  
  prd_df2<-pbmclapply(1:length(prd_df),func_pca,mc.cores=100)
  names(prd_df2)<-names(prd_df)
  
  
  # for(p in 1:length(prd_df2)){
  #   if(ncol(prd_df2[[p]])>0.5*ncol(prd_df[[p]])){prd_df[[p]]=prd_df2[[p]]}
  # }
  prd_df=prd_df2
  prd_df_pca=prd_df2
  
  ####################################
  ###?????????????????????
  prd_mean<-list()
  for(i in 1:length(prd_df)){
    prd_mean[[i]]<-as.data.frame(prd_df[[i]][,1])
    row.names(prd_mean[[i]])<-row.names(prd_df[[i]])
    for(j in 1:length(prd_df[[i]][,1])){
      a1<-sort(as.numeric(prd_df[[i]][j,]))
      prd_mean[[i]][j,1]<-mean(as.numeric(a1))
      #prd_mean[[i]][j,1]<-mean(as.numeric(a1[,ceiling(length(a1)/5):floor(length(a1)/5*4)]))
    }
    names(prd_mean)[i]<-names(prd_df)[i]
    # prd_mean[[i]]<-prd_mean[[i]]/sum(prd_mean[[i]])
  }
  
  ##????????????????????????
  prd_median<-list()
  for(i in 1:length(prd_df)){
    prd_median[[i]]<-as.data.frame(prd_df[[i]][,1])
    row.names(prd_median[[i]])<-row.names(prd_df[[i]])
    for(j in 1:length(prd_df[[i]][,1])){
      a1<-sort(as.numeric(prd_df[[i]][j,]))
      prd_median[[i]][j,1]<-median(as.numeric(a1))
      #prd_mean[[i]][j,1]<-mean(as.numeric(a1[,ceiling(length(a1)/5):floor(length(a1)/5*4)]))
    }
    names(prd_median)[i]<-names(prd_df)[i]
    # prd_mean[[i]]<-prd_mean[[i]]/sum(prd_mean[[i]])
  }
  
  
  qua_sd<-list()
  qua_mean<-list()
  qua_noNA.list<-list()
  
  func_bagging<-function(j){
    print(paste('????????????',j,'???'))
    #df???bagging????????????????????????
    mean1<-prd_mean[[j]]
    mean2<-prd_mean[[j]]
    ##??????????????????
    sd2<-prd_mean[[j]]
    df1<-prd_df[[j]]
    qua_noNA<-list()
    for(i in 1:length(df1[,1])){
      ##Step2
      # df1[is.na(df1)]=0
      # noNA=df1[i,]
      noNA<-df1[i,which(!(is.na(df1[i,])))]
      
      
      noNA_sort<-sort(as.numeric(noNA))
      
      noNA=rbind(noNA,noNA)
      noNA<-noNA[,intersect(colnames(prd_df_pca[[j]]),colnames(noNA))]
      
      noNA=as.data.frame(t(noNA))
      
      ##Step3
      ##??????????????????????????????
      mean2[i,1]<-mean(as.numeric(noNA[,1]))
      sd2[i,1]<-sd(as.numeric(noNA[,1]))
      qua_noNA[[i]]<-noNA
      names(qua_noNA)[i]<-row.names(df1)[i]
      ##???mean1?????????mean2????????????????????????
    }
    
    return(list(sd2,mean2,qua_noNA))
  }
  
  all(names(prd_mean)==names(prd_df_pca))
  
  results_debagging<-pbmclapply(1:length(prd_mean),func_bagging,mc.cores=40)
  names(results_debagging)<-names(prd_mean)
  
  for (j in 1:length(results_debagging)) {
    qua_sd[[j]]<-results_debagging[[j]][[1]]
    qua_mean[[j]]<-results_debagging[[j]][[2]]
    qua_noNA.list[[j]]<-results_debagging[[j]][[3]]
    names(qua_mean)[j]<-names(results_debagging)[j]
    names(qua_noNA.list)[j]<-names(results_debagging)[j]
    names(qua_sd)[j]<-names(results_debagging)[j]
  }
  
  
  ##Step6
  ##????????????
  qua_noNA.listFB<-qua_noNA.list
  qua_meanFB<-qua_mean
  #qua_noNA.list<-qua_noNA.listFB
  for(j in 1:length(qua_noNA.listFB)){
    qua_meanFB[[j]][,1]<-qua_mean[[j]][,1]/sum(qua_mean[[j]][,1])
  }
  
  qua_meanFB_df<-qua_meanFB[[1]]
  colnames(qua_meanFB_df)[1]<-names(qua_meanFB)[1]
  
  for(j in 2:length(qua_meanFB)){
    
    qua_meanFB[[j]]<-as.data.frame(qua_meanFB[[j]])
    qua_meanFB_df<-merge(qua_meanFB_df, qua_meanFB[[j]], by = "row.names", all = TRUE)
    row.names(qua_meanFB_df)<-qua_meanFB_df[,1]
    qua_meanFB_df<-qua_meanFB_df[,-1]
    colnames(qua_meanFB_df)[j]<-names(qua_meanFB)[j]
  }
  qua_meanFB_df2<-qua_meanFB_df
  qua_meanFB_df<-as.data.frame(lapply(qua_meanFB_df2,as.numeric))
  row.names(qua_meanFB_df)<-row.names(qua_meanFB_df2)
  names(qua_meanFB_df)<-names(qua_meanFB)
  
  qua_meanFB_dfFB<-qua_meanFB_df
  for (j in 1:length(qua_meanFB_df[1,])) {
    qua_meanFB_df[,j][is.na(qua_meanFB_df[,j])]<-0}
  
  
  prd_after<-qua_meanFB_df[,sort(names(qua_meanFB_df))]
  
  return(prd_after)
}






cosine_screen_LowToHigh<-function(S,B){
  B=as.data.frame(B)
  B=apply(B,1,sum)
  set.seed(123)
  # B=B[B>0]
  S=S[intersect(names(B),row.names(S)),]
  B=B[intersect(names(B),row.names(S))]
  # 
  ##????????????second_FC??????????????????????????????????????????marker
  second_FC=c()
  for (i in 1:nrow(S)) {
    second_FC=c(second_FC,(sort(S[i,],decreasing = TRUE)[1]+0.01)/(sort(S[i,],decreasing = TRUE)[2]+0.01))
  }
  
  
  S2=as.data.frame(S)
  S2[,'secondFC_name']=names(second_FC)
  S2[,'secondFC_value']=second_FC
  
  unique_name=unique(names(second_FC))
  
  marker_screen=c()
  cos_list=list()
  for (n in 1:length(unique_name)) {
    
    ##???????????????marker???????????????
    second_FC_df_in=S2[S2[,'secondFC_name'] %in% unique_name[n],]
    second_FC_df_in[,'bulk_exp']=B[row.names(second_FC_df_in)]
    second_FC_df_in[,'gene_name']=row.names(second_FC_df_in)
    
    # second_FC_df_in=as.data.frame(second_FC_df_in)
    ##???secondFC??????????????????
    second_FC_df_in=arrange(second_FC_df_in,dplyr::desc(secondFC_value))
    row.names(second_FC_df_in)=second_FC_df_in[,'gene_name']
    
    
    ##?????????marker?????????????????????exp???
    mark_exp=second_FC_df_in[,c(unique_name[n],'bulk_exp')]
    mark_exp=mark_exp[mark_exp[,2]>0,]
    
    a=0
    b=0
    # cor1=cor(mark_exp[,1],mark_exp[,2],method='spearman')
    cor1=lsa::cosine(mark_exp[,1],mark_exp[,2])
    for (v in nrow(mark_exp):1) {
      if(nrow(mark_exp)<10){break}
      mark_exp_in=mark_exp[-v,]
      # cor_in=cor(mark_exp_in[,1],mark_exp_in[,2],method='spearman')
      cor_in=lsa::cosine(mark_exp_in[,1],mark_exp_in[,2])
      
      if(cor_in*nrow(mark_exp_in)>cor1*nrow(mark_exp)){
        mark_exp=mark_exp[-v,]
        cor1=cor_in
        # a=a+1
      }#else{b=b+1}
    }
    cos_list[[n]]=cor1
    marker_screen=c(marker_screen,row.names(mark_exp))
  }
  names(cos_list)=unique_name
  
  S_in=S[marker_screen,]
  data=as.data.frame(S_in)
  # data[,'query']=B[marker_screen]
  B_in=B[row.names(S_in)]
  
  list_cos=cos_list
  list_batch=list(as.matrix(S_in),B_in)
  return(list(list_cos,list_batch))
}


