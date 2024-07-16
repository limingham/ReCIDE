

library(parallel)
library(pbmcapply)
source("~/SWORD/其他模型benchmark测试/CIBERSORT_results/Cibersort_com.R")
####PFC
EXP<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/EXP_and_KEY/EXP.rds")
com_ref<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/ref_data/ref_list.FC_screen.rds")


EXP=EXP[substr(names(com_ref),1,nchar(names(com_ref))-8)]

all(names(EXP)==substr(names(com_ref),1,nchar(names(com_ref))-8))
##com_prd
for(i in 1:length(EXP)){EXP[[i]]<-cbind(EXP[[i]],EXP[[i]])}

fun<-function(x){
  exp<-EXP[[x]]
  ref<-com_ref[[x]]
  results<-CIBERSORT(ref,exp,perm = 100,QN=TRUE)
  return(results)
}

com_results <- pbmclapply(1:length(EXP), fun, mc.cores = 50)
names(com_results)<-names(EXP)
saveRDS(com_results,file='~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_com/CIBEROSRT_com_output.rds')





cibersort_com_result <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_com/CIBEROSRT_com_output.rds")
##comprd结果整合
for (i in 1:length(cibersort_com_result)) {
  cibersort_com_result[[i]]<-as.data.frame(cibersort_com_result[[i]][1,])
  colnames(cibersort_com_result[[i]])<-names(cibersort_com_result)[i]
}
prd<-cibersort_com_result


df_merge <- prd[[1]]
for(j in 2:length(prd)){
  df_merge<-merge(df_merge, prd[[j]], by = "row.names", all = TRUE)
  row.names(df_merge)<-df_merge[,1]
  df_merge<-df_merge[,-1]
  colnames(df_merge)[j]<-names(prd)[j]
}

df_merge<-subset(df_merge,row.names(df_merge)!= 'RMSE')
df_merge<-subset(df_merge,row.names(df_merge)!= 'P-value')
df_merge<-subset(df_merge,row.names(df_merge)!= 'Correlation')
df_merge<-df_merge[,sort(colnames(df_merge))]

com_cibersort<-as.data.frame(lapply(df_merge,as.numeric))
row.names(com_cibersort)<-row.names(df_merge)



prd_com<-com_cibersort



saveRDS(prd_com,file = '~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_com/prd_com_df.rds')