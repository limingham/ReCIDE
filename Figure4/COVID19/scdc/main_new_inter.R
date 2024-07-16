###################################TNBC26
library(Seurat)
library(fastSave)
library(pbmcapply)
library(SCDC)
ref_data<-readRDS.lbzip2('~/scRNA_Seq_data/PBMC/Stephenson_seurat.rdsFS',n.cores = 200)

bulk.mtx<- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/counts_data_query.rds")
rn=intersect(row.names(ref_data),row.names(bulk.mtx))
ref_data=ref_data[rn,]


fdata_bulk <- rownames(bulk.mtx)
pdata_bulk <- cbind(cellname = colnames(bulk.mtx), subjects = colnames(bulk.mtx))
eset_bulk <- getESET(bulk.mtx, fdata = fdata_bulk, pdata = pdata_bulk)



##orig.ident在不同情况下需要改

set.seed(123456)
ref_data_sum=apply(ref_data@assays[["RNA"]]@counts,2,sum)
ref_data_name=names(ref_data_sum[ref_data_sum>10])

ref_data=ref_data[,ref_data_name]
if(ncol(ref_data)>65000){
sample_data=sample(colnames(ref_data), 60000, replace = FALSE, prob = NULL)
ref_data=ref_data[,sample_data]}

data_sc=as.data.frame(as.matrix(ref_data@assays[["RNA"]]@counts))

# source("~/SWORD/SCDC/SCDC-master/R/ENSEMBLE.R")

fdata_sc <- rownames(data_sc)
pdata_sc <- cbind(cellname = colnames(data_sc), 
                  subjects = as.character(ref_data@meta.data[["patient_id"]]),
                  celltype_major=as.character(ref_data@meta.data[["true"]]))

eset_sc <- getESET(data_sc, fdata = fdata_sc, pdata = pdata_sc)
# eset_sc@phenoData@data[,'celltype_major']=ref_data@meta.data[["cg_cov"]]

eset_sc_qc <- SCDC_qc(eset_sc, ct.varname = "celltype_major", sample = "subjects", scsetname = "PBMC",
                      ct.sub = names(table(eset_sc@phenoData@data[["celltype_major"]])), qcthreshold = 0.7) 

# 
prop_results <- SCDC_prop(eset_bulk, eset_sc_qc[["sc.eset.qc"]], ct.varname = "celltype_major", sample = "subjects",
                          ct.sub = unique(ref_data@meta.data[["true"]]))

saveRDS(prop_results,file = '~/ReCIDE/应用_前二_新_inter/COVID19/scdc/SCDC_output_60000.rds')

prd_df=as.data.frame(t(prop_results[["prop.est.mvw"]]))


saveRDS(prd_df,file='~/ReCIDE/应用_前二_新_inter/COVID19/scdc/SCDC_prd_df_60000.rds')
