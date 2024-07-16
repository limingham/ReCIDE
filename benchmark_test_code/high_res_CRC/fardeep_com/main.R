library(FARDEEP)
library(pbmcapply)

####comPBMC
com_ref <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/ref_data/ref_com.FC_screen.rds")
com_ref=com_ref[sort(names(com_ref))]

EXP <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/EXP_and_KEY/EXP.rds")
EXP=EXP[substr(names(com_ref),1,nchar(names(com_ref))-8)]

all(names(EXP)==substr(names(com_ref),1,nchar(names(com_ref))-8))


##com_prd
func_FARDEEP<-function(i){
  EXP_in<-cbind(EXP[[i]],EXP[[i]])
  # REF<-ref_list[[names(EXP)[i]]]
  REF<-com_ref[[i]]
  FD_results_in<-fardeep(REF, EXP_in)
  return(FD_results_in)
}

FARDEEP_results<-pbmclapply(1:length(EXP),func_FARDEEP,mc.cores = 100)
names(FARDEEP_results)<-names(EXP)

saveRDS(FARDEEP_results,file = '~/ReCIDE/benchmark测试_最终/high_res_CRC/fardeep_com/FARDEEP_com_output.rds')

