library(FARDEEP)
library(pbmcapply)

####comPBMC
com_ref<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/ref_data/ref_com.FC_screen.rds")

EXP<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/EXP_and_KEY/EXP.rds")



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

saveRDS(FARDEEP_results,file = '~/ReCIDE/benchmark测试_最终/mult_reference_ER新/fardeep_com/FARDEEP_com_output.rds')

