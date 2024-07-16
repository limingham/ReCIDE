
library(fastSave)
library(DWLS)
library(pbmcapply)

####内并行

EXP<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/EXP_and_KEY/EXP.rds")
EXP=EXP[sort(names(EXP))]
Sig_list <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/DWLS_sep/Sig_ref_list_all.rds")
Sig_list=Sig_list[sort(names(Sig_list))]

source('~/ReCIDE/ReCIDE_main/cosine_screen.R')
# EXP=EXP[unique(seurat_query_healthy@meta.data[["subject"]])]

sep_solDWLS<-list()
fun_DWLS_in<-function(i){
  bulk<-EXP[[i]]
  sep_ref.list<-Sig_list[[i]]
  sep_ref.list[[names(EXP)[i]]]<-NULL
  # sep_ref.list<-ALLPAN_Signature[[names(EXP)[i]]]
  # sep_ref.list=sample(sep_ref.list, 110, replace = FALSE, prob = NULL)
  # for (k in length(sep_ref.list):1) {
  #   if(ncol(sep_ref.list[[k]])<5){sep_ref.list[[k]]<-NULL}
  # }
  # for (k in length(sep_ref.list):1) {
  #   if(length(sep_ref.list[[k]][,1])<length(sep_ref.list[[k]][1,])){sep_ref.list[[k]]<-NULL}
  # }
  
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
sep_solDWLS=pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores=50)
names(sep_solDWLS)=names(EXP)

saveRDS(sep_solDWLS,file='~/ReCIDE/benchmark测试_最终/deone_PAN/DWLS_sep/DWLS_sep_output_LowToHigh.rds')

for (j in 1:length(sep_solDWLS)) {
  sep_solDWLS[[j]]=sep_solDWLS[[j]][[1]]
}
source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
prd_after=ReCIDE_PCA(sep_solDWLS)

prd_after=prd_after[sort(row.names(prd_after)),sort(colnames(prd_after))]

saveRDS(prd_after,file='~/ReCIDE/benchmark测试_最终/deone_PAN/DWLS_sep/prd_sep_LowToHigh_PCA.rds')


