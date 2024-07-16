###################################TNBC26
library(Seurat)
library(fastSave)
library(pbmcapply)
library(SCDC)
ref_data=readRDS.lbzip2('~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/ref_data/ref_all.rdsFS',n.cores=200)
EXP<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/EXP_and_KEY/EXP_combine.rds")
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

rn=intersect(row.names(ref_data),row.names(bulk.mtx))
ref_data=ref_data[rn,]


set.seed(100)

if(ncol(ref_data)>60000){
sample_data=sample(colnames(ref_data), 60000, replace = FALSE, prob = NULL)
ref_data=ref_data[,sample_data]}

ref_data_sum=apply(ref_data@assays[["RNA"]]@counts,2,sum)
ref_data_name=names(ref_data_sum[ref_data_sum>10])

ref_data=ref_data[,ref_data_name]

data_sc=as.data.frame(as.matrix(ref_data@assays[["RNA"]]@counts))

# source("~/SWORD/SCDC/SCDC-master/R/ENSEMBLE.R")

fdata_sc <- rownames(data_sc)
pdata_sc <- cbind(cellname = colnames(data_sc), 
                  subjects = as.character(ref_data@meta.data[["ind_cov"]]),
                  celltype_major=as.character(ref_data@meta.data[["cg_cov"]]))

eset_sc <- getESET(data_sc, fdata = fdata_sc, pdata = pdata_sc)
# eset_sc@phenoData@data[,'celltype_major']=ref_data@meta.data[["cg_cov"]]

eset_sc_qc <- SCDC_qc(eset_sc, ct.varname = "celltype_major", sample = "subjects", scsetname = "PBMC",
                      ct.sub = names(table(eset_sc@phenoData@data[["celltype_major"]])), qcthreshold = 0.7) 

# 
prop_results <- SCDC_prop(eset_bulk, eset_sc_qc[["sc.eset.qc"]], ct.varname = "celltype_major", sample = "subjects",
                          ct.sub = unique(ref_data@meta.data[["cg_cov"]]))

saveRDS(prop_results,file = '~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/scdc/SCDC_output.rds')

prd_df=as.data.frame(t(prop_results[["prop.est.mvw"]]))

key_df <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/EXP_and_KEY/key_combine.rds")
prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]


RMSE_vec=c()
for(i in 1:length(key_df)){
  RMSE_vec[i]=ModelMetrics::rmse(key_df[,i],prd_df[,i])
  # RMSE_vec[i]=cor(key_df[,i],prd_df[,i])
  
}
mean(RMSE_vec)

saveRDS(prd_df,file='~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/scdc/SCDC_prd_df.rds')
