


source("~/ReCIDE/ReCIDE_main/先PCA再层次聚类.R")
source('~/ReCIDE/ReCIDE_main/cosine_screen.R')
####sepPFC
sep_ref.list <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/ref_data/ER_sep_sig.rds")
EXP<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/EXP_and_KEY/EXP.rds")
source("~/SWORD/其他模型benchmark测试/CIBERSORT_results/Cibersort.R")


sep_solFAR<-list()

for (i in 1:length(EXP)){
  
  bulk<-EXP[[i]]
  ref.list_in<-sep_ref.list
  # sep_ref.list[[names(EXP)[i]]]=0
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
  
  sep_solFAR[[i]]<-pbmclapply(1:length(ref.list_in),fun_FAR_in,mc.cores = 100)
  names(sep_solFAR[[i]])=names(ref.list_in)
  
}
names(sep_solFAR)<-names(EXP)

saveRDS(sep_solFAR,file='~/ReCIDE/benchmark测试_最终/mult_reference_ER新/cibersort_sep/sep_CIBERSORT_output_LowToHigh.rds')


