
source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
source("~/ReCIDE/ReCIDE_main/cosine_screen.R")

library(fastSave)
library(DWLS)

Sig_list <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/DWLS_sep/Sig_ref_list_all.rds")
Sig_list=Sig_list[sort(names(Sig_list))]

EXP<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/EXP_and_KEY/EXP.rds")

EXP=EXP[sort(names(EXP))]


sep_solDWLS<-list()
fun_DWLS_in<-function(i){
  bulk<-EXP[[i]]
  sep_ref.list<-Sig_list[[i]]
  sep_ref.list[[names(EXP)[i]]]<-NULL
  
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


sep_solDWLS<-pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores = 110)
names(sep_solDWLS)<-names(EXP)
# }

saveRDS(sep_solDWLS,file='~/ReCIDE/benchmark测试_最终/high_res_CRC/DWLS_sep/DWLS_sep_output.rds')

sep_solDWLS=readRDS('~/ReCIDE/benchmark测试_最终/high_res_CRC/DWLS_sep/DWLS_sep_output.rds')

source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")
source("~/ReCIDE/ReCIDE_main/cosine_screen.R")

# for (i in 1:length(sep_solDWLS)) {
#   for (j in 1:length(sep_solDWLS[[i]])) {
#   sep_solDWLS[[i]][[j]]=sep_solDWLS[[i]][[j]][[1]]
#   }
#   sep_solDWLS[[i]]=sep_solDWLS[[i]][[1]]
# }

# 
for (i in 1:length(sep_solDWLS)) {
  sep_solDWLS[[i]]=sep_solDWLS[[i]][[1]]

}

prd_after=ReCIDE_PCA(sep_solDWLS)


prd_after=prd_after[sort(row.names(prd_after)),sort(colnames(prd_after))]


saveRDS(prd_after,file='~/ReCIDE/benchmark测试_最终/high_res_CRC/DWLS_sep/prd_sep_LowToHigh_PCA.rds')

