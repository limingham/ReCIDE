
library(fastSave)
library(DWLS)


deconvolution_function <- function(Sig_list,EXP_df,n_cores=5,n_celltype=0){
  
  EXP=list()
  for (i in 1:ncol(EXP_df)) {
    EXP[[i]]=as.numeric(EXP_df[,i])
    names(EXP[[i]])=row.names(EXP_df)
  }
  names(EXP)=colnames(EXP_df)
  
  
  sep_solDWLS<-list()
  fun_DWLS_in<-function(i){
    bulk<-EXP[[i]]
    sep_ref.list<-Sig_list
    
    sep_solDWLS_inside<-list()
    cos_list_inside=list()
    for (j in 1:length(sep_ref.list)) {
      ##dwls??????dataframe?????????matrix
      ref1<-as.matrix(sep_ref.list[[j]])
      
      # batch_output<-cosine_screen_HighToLow(ref1,bulk)
      batch_output<-cosine_screen_LowToHigh(ref1,bulk)
      # batch_output<-cosine_screen_GA(ref1,bulk)
      
      tr <- batch_output[[2]]
      # tr<-trimData(ref1,bulk)
      sep_solDWLS_inside[[j]]<-try(solveDampenedWLS(tr[[1]],tr[[2]]), TRUE)
      
      cos_list_inside[[j]]<-batch_output[[1]]
      
      names(sep_solDWLS_inside)[j]<-names(sep_ref.list)[j]
      names(cos_list_inside)[j]<-names(sep_ref.list)[j]
      
      DWLS_output=list(sep_solDWLS_inside,cos_list_inside)
      names(DWLS_output)=c('sep_solDWLS','cos_list')
    }
    return(DWLS_output)
  }
  sep_solDWLS=pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores=n_cores)
  names(sep_solDWLS)=names(EXP)
  
  
  for (j in 1:length(sep_solDWLS)) {
    sep_solDWLS[[j]]=sep_solDWLS[[j]][[1]]
    for(k in length(sep_solDWLS[[j]]):1){
      if(length(sep_solDWLS[[j]][[k]])<n_celltype){sep_solDWLS[[j]][[k]]<-NULL}
    }
  }
  # source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
  prd_df=ReCIDE_PCA(sep_solDWLS)
  
  prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
  
  
  deconvolution_output=list(sep_solDWLS,prd_df)
  names(deconvolution_output)=c('results_list','results_final_df')
  return(deconvolution_output)
}

