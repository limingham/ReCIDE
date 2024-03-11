suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(ggbiplot))
suppressMessages(library(pbmcapply))
suppressMessages(library(mclust))



ReCIDE_PCA=function(sep_solDWLS,method='DWLS'){
  
  for (i in length(sep_solDWLS):1) {
    for (j in length(sep_solDWLS[[i]]):1) {
      if(class(sep_solDWLS[[i]][[j]])[1]=="try-error"){sep_solDWLS[[i]][[j]]<-NULL}
    }
  }
  
  if(method=='FARDEEP'){
    for (i in 1:length(sep_solDWLS)) {
      for (j in 1:length(sep_solDWLS[[i]])) {
        sep_solDWLS[[i]][[j]]<-sep_solDWLS[[i]][[j]][["relative.beta"]][1,]
      }}
    prd<-sep_solDWLS
  }
  ##DWLS结果整合
  # # prd<-sep_solOLS
  if(method=='DWLS'){
    for (i in length(sep_solDWLS):1) {
      for (j in length(sep_solDWLS[[i]]):1) {
        if(class(sep_solDWLS[[i]][[j]])=="try-error"){sep_solDWLS[[i]][[j]]<-NULL}
      }
    }
    prd<-sep_solDWLS
  }
  
  # ##delia整合
  # for (i in 1:length(delia_sep_results)) {
  #   for (j in 1:length(delia_sep_results[[i]])) {delia_sep_results[[i]][[j]]<-delia_sep_results[[i]][[j]][["out"]][,1]
  #   }
  # }
  # prd<-delia_sep_results
  
  
  # cibersort结果整合
  if(method=='CIBERSORT'){
    for (i in 1:length(sep_solDWLS)) {
      for (j in 1:length(sep_solDWLS[[i]])) {
        data_in=sep_solDWLS[[i]][[j]]
        data_in<-as.data.frame(data_in[1,-((ncol(data_in)-2):ncol(data_in))])
        data_in=as.numeric(data_in)
        names(data_in)=colnames(sep_solDWLS[[i]][[j]])[1:(ncol(sep_solDWLS[[i]][[j]])-3)]
        sep_solDWLS[[i]][[j]]=data_in
      }
    }
    prd<-sep_solDWLS
  }
  
  #prd是列表，变为df格式
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
  
  # for (j in 1:length(prd_df)) {
  #   prd_df[[j]][is.na(prd_df[[j]])]<-0
  # }
  # 
  
  #prd_df<-prd_dfFB
  prd_dfFB<-prd_df
  for (j in 1:length(prd_df)) {
    for(i in 1:length(prd_df[[j]][,1])){
      ###0.75的阈值合适吗？
      if(table(is.na(prd_df[[j]][i,]))['FALSE']<0.75*ncol(prd_df[[j]])){
        prd_df[[j]][i,][is.na(prd_df[[j]][i,])]<-0}
    }
  }
  ####如果做PCA的话，从此处导入PCA模块
  ##PCA
  ####PCA模块 
  
  #################
  # for (n in 1:length(prd_df)) 
  func_pca<- function(n){
    
    data_train <- t(prd_df[[n]])
    
    data_train[is.na(data_train)]=0
    
    for (i in ncol(data_train):1) {
      #如果某列中的0值大于行长度的一半，就去掉
      if((nrow(as.data.frame(data_train[data_train[,i]<=0,]))>0.5*nrow(data_train))){
        data_train=data_train[,-i]}
    }
    data_train=as.data.frame(data_train)
    if(ncol(data_train)>=3){
      library(gmodels)
      #将带有训练数据的数据框更改为矩阵
      PCA.OUT=fast.prcomp(data_train, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL)
      
      imp=PCA.OUT[["sdev"]]/sum(PCA.OUT[["sdev"]])
      
      # im=0
      # for (l in 1:length(imp)) {
      #   im=imp[l]+im
      #   if((imp[l]<0.1) | (im>0.9)){
      #     pc_sure=l
      #     break
      #   }
      # }
      
      # if(length(imp)<=10){
        pc_sure=2
      # }else if(length(imp)<=20){
      #   pc_sure=3
      # }else{
      #   pc_sure=4
      # }
      
      pca_n=PCA.OUT$x[,1:pc_sure]
      
      pca_n[,1]=scale(pca_n[,1])
      pca_n[,2]=scale(pca_n[,2])
      
      ##pc百分比小于10%
      ##总贡献大于70%
      ##总PC数小于10
      ##聚类方式1，mclust
      #####################################
      #####################################
      model <- Mclust(pca_n)

      
      # & (model_class[1,2]/nrow(pca_n)>0.5)
      ###想不跳过PCA需要修改此处
      if(model$G !=1){
        model_class=as.data.frame(table(model$classification))
        model_class=model_class[order(model_class[,'Freq'],decreasing = TRUE),]

        if(model_class[1,2]>5){
        if(model_class[1,2]==model_class[2,2]){
          ##如果出现两类聚类数相同
          
          model_class=model_class[model_class[,2]==model_class[1,2],]
          
          use_ind=model$classification[model$classification %in% model_class[,1]]
          data_train=data_train[names(use_ind),]
          
          data_train_list=split(data_train,use_ind) 

          kappa_value=c()
          for (ka in 1:length(data_train_list)) {
            data_in=as.data.frame(t(data_train_list[[ka]]))
            kappa_value[ka]=kappa(data_in)
          }
          
          # kappa<100则认为有共线性程度很小。
          # 100<=kappa<=1000则认为存在中等程度或者较强的共线性。
          # 若kappa>1000则认为存在很严重的共线性

          patient_names=row.names(data_train_list[[which.max(kappa_value)]])
        }else{
          ##如果出现两类聚类数不同
          patient_names=names(model$classification[model$classification==model_class[1,1]])
        }
        }else{
            ##这个else对应的是大于5
          patient_names=colnames(prd_df[[n]])
          }
        
      }else{
        ##如果聚类数仅为1，则不筛选
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
  ###计算预测的均值
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
  
  ##计算预测的中位数
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
    print(paste('现在是第',j,'组'))
    #df是bagging前的预测结果矩阵
    mean1<-prd_mean[[j]]
    mean2<-prd_mean[[j]]
    ##值是随便赋的
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
      ##重新求出取点后的均值
      mean2[i,1]<-mean(as.numeric(noNA[,1]))
      sd2[i,1]<-sd(as.numeric(noNA[,1]))
      qua_noNA[[i]]<-noNA
      names(qua_noNA)[i]<-row.names(df1)[i]
      ##将mean1重铸为mean2后，进行后续操作
    }
    #}
    # qua_sd[[j]]<-sd2
    # qua_mean[[j]]<-mean2
    # qua_noNA.list[[j]]<-qua_noNA
    # names(qua_mean)[j]<-names(prd_mean)[j]
    # names(qua_noNA.list)[j]<-names(prd_mean)[j]
    # names(qua_sd)[j]<-names(prd_mean)[j]
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
  ##风险模型
  qua_noNA.listFB<-qua_noNA.list
  qua_meanFB<-qua_mean
  #qua_noNA.list<-qua_noNA.listFB
  for(j in 1:length(qua_noNA.listFB)){
    # index_minus<-matrix(nrow = length(qua_noNA.listFB[[j]]), ncol = 1)
    # index_minus<-as.data.frame(index_minus)
    # for(i in 1:length(qua_noNA.listFB[[j]])){
    #   qua_noNA.listFB[[j]][[i]]<-as.data.frame(qua_noNA.listFB[[j]][[i]])
    #   ##此处qua_noNA.list和qua_mean是否是真正对应？
    #   qua_noNA.listFB[[j]][[i]][(length(qua_noNA.listFB[[j]][[i]][1,])+1),1]<-qua_mean[[j]][i,1]
    #   ##排序并筛选
    #   sortdata<-sort(qua_noNA.listFB[[j]][[i]][,1])
    #   index<-which(sortdata==qua_mean[[j]][i,1])
    #   #index<-floor(mean(index))
    #   if((index < 0.7*length(sortdata)) & (index > 0.3*length(sortdata))){index_minus[i,1]<-1
    #   }else{index_minus[i,1]<-1+mean(min(abs(0.3*length(sortdata)-index),abs(0.7*length(sortdata)-index)))/length(sortdata)
    # }}
    # 
    # x=1/sum(qua_mean[[j]][,1] * index_minus)
    # qua_meanFB[[j]][,1]<-qua_mean[[j]][,1] * index_minus * x
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
  
  # qua_meanFB_df<-as.data.frame(qua_meanFB)
  # colnames(qua_meanFB_df)<-names(qua_meanFB)
  prd_after<-qua_meanFB_df[,sort(names(qua_meanFB_df))]
  
  return(prd_after)
}
# save(prd_after,file='~/SWORD/新想法/贝叶斯回归/prd_after_CRC.RData')

