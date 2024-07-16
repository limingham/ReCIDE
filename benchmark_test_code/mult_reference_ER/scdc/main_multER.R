library(Seurat)
library(fastSave)
library(pbmcapply)
library(SCDC)
ER.dgedata_150=readRDS.lbzip2('~/ReCIDE/benchmark测试_最终/mult_reference_ER新/ref_data/ref_all.rdsFS',n.cores = 300)



#######################构建bulk
# TNBC.dgedata_150.list<-SplitObject(TNBC.dgedata_150,split.by = 'orig.ident')

EXP<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/EXP_and_KEY/EXP.rds")
bulk.mtx<-as.data.frame(EXP[[1]])
for (i in 2:length(EXP)) {
  bulk.mtx2<-as.data.frame(EXP[[i]])
  bulk.mtx<-cbind(bulk.mtx,bulk.mtx2)
  
}
bulk.mtx<-as.matrix(bulk.mtx)
colnames(bulk.mtx)<-names(EXP)

gene_exp=apply(bulk.mtx,1,sum)
bulk.mtx=bulk.mtx[gene_exp>0,]

fdata_bulk <- rownames(bulk.mtx)
pdata_bulk <- cbind(cellname = colnames(bulk.mtx), subjects = colnames(bulk.mtx))
eset_bulk <- getESET(bulk.mtx, fdata = fdata_bulk, pdata = pdata_bulk)





##orig.ident在不同情况下需要改

ref_data_in_list=SplitObject(ER.dgedata_150,split.by = "research")

# eset_sc_qc.list=list()
# for (k in 1:length(ref_data_in_list))
QCK<-function(k){
  set.seed(100)
  ref_data_in=ref_data_in_list[[k]]
  
  rn=intersect(row.names(ref_data_in),row.names(bulk.mtx))
  ref_data_in=ref_data_in[rn,]
  
  
  if(ncol(ref_data_in)>60000){
    sample_data=sample(colnames(ref_data_in), 60000, replace = FALSE, prob = NULL)
    ref_data_in=ref_data_in[,sample_data]}
  # ref_data_in <- sample_seob(ref_data_in,sp.total =65000,group.by='orig.ident')
  
  ref_data_sum=apply(ref_data_in@assays[["RNA"]]@counts,2,sum)
  ref_data_name=names(ref_data_sum[ref_data_sum>10])
  
  ref_data_in=ref_data_in[,ref_data_name]
  
  data_sc=as.data.frame(as.matrix(ref_data_in@assays[["RNA"]]@counts))
  
  
  fdata_sc <- rownames(data_sc)
  pdata_sc <- cbind(cellname = colnames(data_sc), 
                    subjects = ref_data_in@meta.data[["batch"]],
                    celltype_major=ref_data_in@meta.data[["usetype"]])
  
  eset_sc <- getESET(data_sc, fdata = fdata_sc, pdata = pdata_sc)
  # eset_sc@phenoData@data[,'celltype_major']=ref_data@meta.data[["celltype_major"]]
  
  eset_sc_qc <- SCDC_qc(eset_sc, ct.varname = "celltype_major", sample = "subjects", scsetname = "TNBC",
                        ct.sub = names(table(eset_sc@phenoData@data[["celltype_major"]])), qcthreshold = 0.7) 
  
  # eset_sc_qc.list[[k]]=eset_sc_qc$sc.eset.qc
  # names(eset_sc_qc.list)[k]=names(ref_data_in_list)[k]
  return(eset_sc_qc$sc.eset.qc)
}

eset_sc_qc.list=pbmclapply(1:length(ref_data_in_list),QCK,mc.cores=2)

names(eset_sc_qc.list)=names(ref_data_in_list)

saveRDS(eset_sc_qc.list,file = '~/ReCIDE/benchmark测试_最终/mult_reference_ER新/scdc/SCDC_qc_results.rds')


source("~/SWORD/SCDC/SCDC-master/R/ENSEMBLE.R")

# This might take several minutes, depending on the search.length set by user.
prop_results <- SCDC_ENSEMBLE2(bulk.eset = eset_bulk, sc.eset.list = eset_sc_qc.list, ct.varname = "celltype_major",
                               sample = "subjects", truep = NULL, 
                               ct.sub =  unique(ER.dgedata_150@meta.data[,'usetype']), search.length = 0.01, grid.search = T)

saveRDS(prop_results,file = '~/ReCIDE/benchmark测试_最终/mult_reference_ER新/scdc/SCDC_output.rds')

SCDC_output <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/scdc/SCDC_output.rds")
prop_df=wt_prop(SCDC_output[["w_table"]]['inverse SSE',1:2], SCDC_output$prop.only)
# prop_df=wt_prop(SCDC_output[["w_table"]]['Spearman_Y',1:2], SCDC_output$prop.only)

prop_df=as.data.frame(t(prop_df))
saveRDS(prop_df,file = '~/ReCIDE/benchmark测试_最终/mult_reference_ER新/scdc/SCDC_prop_df.rds')



