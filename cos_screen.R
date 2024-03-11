library(DWLS)
library(dplyr)
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(ggbiplot))
suppressMessages(library(pbmcapply))
suppressMessages(library(mclust))
# 

cosine_screen_LowToHigh<-function(S,B){
  B=as.data.frame(B)
  B=apply(B,1,sum)
  set.seed(123)
  # B=B[B>0]
  S=S[intersect(names(B),row.names(S)),]
  B=B[intersect(names(B),row.names(S))]
  # 
  ##通过计算second_FC判断每个基因是哪种细胞类型的marker
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
    
    ##只留下这种marker基因的某行
    second_FC_df_in=S2[S2[,'secondFC_name'] %in% unique_name[n],]
    second_FC_df_in[,'bulk_exp']=B[row.names(second_FC_df_in)]
    second_FC_df_in[,'gene_name']=row.names(second_FC_df_in)
    
    ##按secondFC的值降序排列
    second_FC_df_in=arrange(second_FC_df_in,desc(secondFC_value))
    row.names(second_FC_df_in)=second_FC_df_in[,'gene_name']
    
    
    ##只留下marker基因对应的列和exp列
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