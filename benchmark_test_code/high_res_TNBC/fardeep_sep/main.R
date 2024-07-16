library(FARDEEP)
library(pbmcapply)

####sepPBMC
load("~/SWORD/单细胞图谱/参考集/整理EXP/EXP_PBMC.RData")
EXP <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/EXP_and_KEY/EXP.rds")
sep_ref.list<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/ref_data/Sig_ref_list_all.rds")
source("~/ReCIDE/ReCIDE_main/先PCA再层次聚类.R")
source('~/ReCIDE/ReCIDE_main/cosine_screen.R')


sep_ref.list=sep_ref.list[sort(names(sep_ref.list))]
EXP=EXP[sort(names(EXP))]

func_FARDEEP<-function(i){
  EXP_in<-cbind(EXP[[i]],EXP[[i]])
  bulk=EXP_in
  
  ref.list_in=sep_ref.list[[i]]
  ref.list_in[[names(EXP)[i]]]<-NULL
  
  FD_results_in<-list()
  for (j in 1:length(ref.list_in)) {
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
    
    FD_results_in[[j]]<-fardeep(ref_df, EXP_in)
  }
  names(FD_results_in)<-names(ref.list_in)
  
  return(FD_results_in)
}

FARDEEP_results<-pbmclapply(1:length(EXP),func_FARDEEP,mc.cores = 110)
names(FARDEEP_results)<-names(EXP)


saveRDS(FARDEEP_results,file='~/ReCIDE/benchmark测试_最终/high_res_TNBC/fardeep_sep/FARDDEP_sep_output_LowToHigh.rds')

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


library(FARDEEP)
####sepPFC
sep_ref.list<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/ref_data/PFC_sep_sig.rds")

EXP <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/EXP_and_key/EXP_PFC.rds")

all(names(EXP)==names(sep_ref.list))

source("~/ReCIDE/ReCIDE_main/先PCA再层次聚类.R")
source('~/ReCIDE/ReCIDE_main/cosine_screen.R')
sep_solFAR<-list()

for (i in 1:length(EXP)) {
  
  bulk<-EXP[[i]]
  ref.list_in<-sep_ref.list
  ref.list_in[[names(EXP)[i]]]<-NULL
  
  FD_results_in<-list()
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
    
    FD_results_in<-fardeep(ref_df, query_df)
    return(FD_results_in)
  }
  
  sep_solFAR[[i]]<-pbmclapply(1:length(ref.list_in),fun_FAR_in,mc.cores = 110)
  names(sep_solFAR[[i]])=names(ref.list_in)
  
}
names(sep_solFAR)<-names(EXP)

saveRDS(sep_solFAR,file='~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/fardeep_sep/FARDDEP_sep_output_LowToHigh.rds')
