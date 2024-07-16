
library(dplyr)

source('~/ReCIDE/ReCIDE_main/cosine_screen.R')
####sepPFC
sep_ref.list <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/ref_data/Sig_ref_list_all.rds")
EXP<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/EXP_and_KEY/EXP.rds")
source("~/SWORD/其他模型benchmark测试/CIBERSORT_results/Cibersort.R")


sep_ref.list=sep_ref.list[sort(names(sep_ref.list))]
EXP=EXP[sort(names(EXP))]



sep_solFAR<-list()

for (i in 1:length(EXP)){
  
  bulk<-EXP[[i]]
  ref.list_in<-sep_ref.list[[i]]
  ref.list_in[[names(EXP)[i]]]<-NULL
  # sep_ref.list<-ALLPAN_Signature[[names(EXP)[i]]]
  # sep_ref.list=sample(sep_ref.list, 110, replace = FALSE, prob = NULL)
  
  ##根据细胞类型数筛选一次
  
  CB_results_in<-list()
  fun_FAR_in<-function(j){
    ##dwls??????dataframe?????????matrix
    ref1<-as.matrix(ref.list_in[[j]])
    
    # batch_output<-cosine_screen_HighToLow(ref1,bulk)
    batch_output<-cosine_screen_LowToHigh(ref1,bulk)
    # batch_output<-cosine_screen_GA(ref1,bulk)
    tr <- batch_output[[2]]
    # tr<-trimData(ref1,bulk)
    query_df=as.data.frame(tr[[2]])
    query_df=cbind(query_df,query_df)
    
    ref_df=as.data.frame(tr[[1]])
    CB_results_in<-CIBERSORT(ref_df,query_df,perm = 100,QN=TRUE)
    
    return(CB_results_in)
  }
  
  sep_solFAR[[i]]<-pbmclapply(1:length(ref.list_in),fun_FAR_in,mc.cores = 50)
  names(sep_solFAR[[i]])=names(ref.list_in)
  
}
names(sep_solFAR)<-names(EXP)

saveRDS(sep_solFAR,file='~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_sep/sep_CIBERSORT_output_LowToHigh.rds')

rm(list=ls())
gc()



sep_solDWLS <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_sep/sep_CIBERSORT_output_LowToHigh.rds")
source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
# source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA_mclust2.R")
prd_after=ReCIDE_PCA(sep_solDWLS,method='CIBERSORT')
saveRDS(prd_after,file='~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_sep/prd_sep_LowToHigh_PCA2.rds')

prd_df<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/cibersort_com/prd_com_df.rds")
KEY <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/EXP_and_KEY/KEY.rds")

prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
KEY=KEY[sort(row.names(KEY)),sort(colnames(KEY))]

rmse_vec=c()
for (i in 1:ncol(KEY)) {
 rmse_vec[i] = ModelMetrics::rmse(as.numeric(KEY[,i]),as.numeric(prd_df[,i]))
  
}

mean(rmse_vec)
