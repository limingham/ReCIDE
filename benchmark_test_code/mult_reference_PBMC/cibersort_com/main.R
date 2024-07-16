

library(parallel)
library(pbmcapply)
source("~/SWORD/其他模型benchmark测试/CIBERSORT_results/Cibersort_com.R")
####PFC
EXP<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/EXP_and_KEY/EXP_PBMC.rds")
com_ref <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/ref_data/ref_com.FC_screen.rds")

# EXP=EXP[substr(names(com_ref),1,nchar(names(com_ref))-8)]

# all(names(EXP)==substr(names(com_ref),1,nchar(names(com_ref))-8))
##com_prd
for(i in 1:length(EXP)){EXP[[i]]<-cbind(EXP[[i]],EXP[[i]])}

fun<-function(x){
  exp<-EXP[[x]]
  ref<-com_ref
  results<-CIBERSORT(ref,exp,perm = 100,QN=TRUE)
  return(results)
}

com_results <- pbmclapply(1:length(EXP), fun, mc.cores = 50)
names(com_results)<-names(EXP)
saveRDS(com_results,file='~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/cibersort_com/CIBEROSRT_com_output.rds')

