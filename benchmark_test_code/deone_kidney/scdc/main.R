###################################TNBC26
library(Seurat)
library(fastSave)
library(pbmcapply)
library(SCDC)
# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
combine.dgedata_150=readRDS.lbzip2('~/ReCIDE/benchmark测试_最终/deone_kidney/ref_data/ref_all.rdsFS',n.cores = 200)
patient_names=unique(combine.dgedata_150@meta.data['patient'])


EXP<- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/EXP_and_KEY/EXP_kidney.rds")
bulk.mtx<-as.data.frame(EXP[[1]])
for (i in 2:length(EXP)) {
  bulk.mtx2<-as.data.frame(EXP[[i]])
  bulk.mtx<-cbind(bulk.mtx,bulk.mtx2)
  
}
bulk.mtx<-as.matrix(bulk.mtx)
colnames(bulk.mtx)<-names(EXP)

fdata_bulk <- rownames(bulk.mtx)
pdata_bulk <- cbind(cellname = colnames(bulk.mtx), subjects = colnames(bulk.mtx))
eset_bulk <- getESET(bulk.mtx, fdata = fdata_bulk, pdata = pdata_bulk)





func_in<-function(f){
  ##orig.ident在不同情况下需要改

  ref_data<-subset(combine.dgedata_150,patient != patient_names[f,1])


  if(ncol(ref_data)>60000){
    sample_data=sample(colnames(ref_data), 60000, replace = FALSE, prob = NULL)
    ref_data=ref_data[,sample_data]}
  
  
  # source("~/SWORD/SCDC/SCDC-master/R/ENSEMBLE.R")
  ref_data_sum=apply(ref_data@assays[["RNA"]]@counts,2,sum)
  ref_data_name=names(ref_data_sum[ref_data_sum>10])
  ref_data=ref_data[,ref_data_name]
  
  data_sc=as.data.frame(as.matrix(ref_data@assays[["RNA"]]@counts))
  
  fdata_sc <- rownames(data_sc)
  pdata_sc <- cbind(cellname = colnames(data_sc), 
                    subjects = as.character(ref_data@meta.data[["patient"]]),
                    celltype_major=as.character(ref_data@meta.data[["subclass.l1"]]))
  
  eset_sc <- getESET(data_sc, fdata = fdata_sc, pdata = pdata_sc)
  # eset_sc@phenoData@data[,'celltype_major']=ref_data@meta.data[["cg_cov"]]
  
  rm(data_sc)
  gc()
  
  eset_sc_qc <- SCDC_qc(eset_sc, ct.varname = "celltype_major", sample = "subjects", scsetname = "PBMC",
                        ct.sub = names(table(eset_sc@phenoData@data[["celltype_major"]])), qcthreshold = 0.7) 
  
  # 
  
  rm(eset_sc)
  gc()
  
  prop_results <- SCDC_prop(eset_bulk, eset_sc_qc[["sc.eset.qc"]], ct.varname = "celltype_major", sample = "subjects",
                            ct.sub = unique(ref_data@meta.data[["subclass.l1"]]))
  
  return(prop_results)
}

prop_results.list=pbmclapply(1:nrow(patient_names),func_in,mc.cores = 4)
names(prop_results.list)=patient_names[,1]

saveRDS(prop_results.list,file = '~/ReCIDE/benchmark测试_最终/deone_kidney/scdc/SCDC_output.rds')



de_1=as.data.frame(t(SCDC_output[[1]][["prop.est.mvw"]])[,names(SCDC_output)[1]])
for(i in 2:length(SCDC_output)){
  de_in=as.data.frame(t(SCDC_output[[i]][["prop.est.mvw"]])[,names(SCDC_output)[i]])
  
  de_1<-merge(de_1, de_in, by = "row.names", all = TRUE)
  row.names(de_1)<-de_1[,1]
  de_1<-de_1[,-1]
}
colnames(de_1)=names(SCDC_output)
# de_1['NEU',]=0
saveRDS(de_1,file='~/ReCIDE/benchmark测试_最终/deone_kidney/scdc/prd_df.rds')

# prd_df=as.data.frame(t(prop_results[["prop.est.mvw"]]))
