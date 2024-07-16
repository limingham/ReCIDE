library(FARDEEP)
library(pbmcapply)

####comPBMC
com_ref <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/ref_data/ref_com.FC_screen.rds")

EXP <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/EXP_and_KEY/EXP_combine.rds")



##com_prd
func_FARDEEP<-function(i){
  EXP_in<-cbind(EXP[[i]],EXP[[i]])
  # REF<-ref_list[[names(EXP)[i]]]
  REF<-com_ref
  FD_results_in<-fardeep(REF, EXP_in)
  return(FD_results_in)
}

FARDEEP_results<-pbmclapply(1:length(EXP),func_FARDEEP,mc.cores = 100)
names(FARDEEP_results)<-names(EXP)

saveRDS(FARDEEP_results,file = '~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/fardeep_com/FARDEEP_com_output.rds')




prd<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/fardeep_com/FARDEEP_com_output.rds")

prd[[1]]<-as.data.frame(prd[[1]][["relative.beta"]])
df_merge<-as.data.frame(t(prd[[1]])[,1])
colnames(df_merge)[1]<-names(prd)[1]

for(j in 2:length(prd)){
  
  prd[[j]]<-as.data.frame(prd[[j]][["relative.beta"]])
  merge_in=as.data.frame(t(prd[[j]])[,1])
  
  df_merge<-merge(df_merge,merge_in, by = "row.names", all = TRUE)
  row.names(df_merge)<-df_merge[,1]
  df_merge<-df_merge[,-1]
  colnames(df_merge)[j]<-names(prd)[j]
}

prd_com<-df_merge[,sort(names(df_merge))]


prd_com=prd_com[sort(row.names(prd_com)),sort(colnames(prd_com))]


saveRDS(prd_com,file = '~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/fardeep_com/prd_com_df.rds')
