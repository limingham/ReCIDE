library(FARDEEP)
library(pbmcapply)

####sepPBMC
Sig_list <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/ref_data/Sig_ref_list_all.rds")
Sig_list=Sig_list[sort(names(Sig_list))]

EXP<- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/EXP_and_KEY/EXP_kidney.rds")
EXP[["28-10051"]]<-NULL
EXP[["33-10006"]]<-NULL
EXP=EXP[sort(names(EXP))]

source("~/ReCIDE/ReCIDE_main/先PCA再层次聚类.R")
source('~/ReCIDE/ReCIDE_main/cosine_screen.R')

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


saveRDS(FARDEEP_results,file='~/ReCIDE/benchmark测试_最终/deone_kidney/fadeep_sep/FARDDEP_sep_output_LowToHigh.rds')

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


sep_solDWLS <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/fadeep_sep/FARDDEP_sep_output_LowToHigh.rds")
source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
# source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA_mclust2.R")
prd_after=ReCIDE_PCA(sep_solDWLS,method='FARDEEP')
saveRDS(prd_after,file='~/ReCIDE/benchmark测试_最终/deone_kidney/fadeep_sep/prd_df.rds')
