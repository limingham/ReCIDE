library(Seurat)
library(fastSave)
library(pbmcapply)
library(SCDC)
combined.data<-readRDS.lbzip2('~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/Seurat_haoANDstep_combine.rdsFS',n.cores = 200)
EXP_combine <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/EXP_and_KEY/EXP_PBMC.rds")

cell_names=intersect(row.names(combined.data),names(EXP_combine[[1]]))
combined.data=combined.data[cell_names,]


patient_id<-as.data.frame(table(combined.data@meta.data[["batch"]]))
patient_id<-patient_id[patient_id[,2]>500,]
combined.data<-subset(combined.data,batch %in% patient_id[,1])

cell_id=as.data.frame(table(combined.data@meta.data[["cell_label"]]))
cell_id=cell_id[cell_id[,2]>20,]
combined.data<-subset(combined.data,cell_label %in% cell_id[,1])


gc()


rowsum<-rowSums(combined.data@assays[["RNA"]]@counts)
rowsum<-as.data.frame(rowsum)
rowsum2<-subset(rowsum,rowsum>0)
combined.data<-combined.data[row.names(rowsum2),]
combined.data <- NormalizeData(combined.data, verbose = FALSE)


Idents(combined.data)<-'cell_label'
COSG_markers <- cosg(
  combined.data,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)

marker_list = list()
for (i in 1:length(COSG_markers[["names"]])) {
  marker_list[[i]]=COSG_markers[["names"]][,i]
}
#list_marker = marker_list

unique_cell_type= colnames(COSG_markers[["names"]])
#cell_type_unique = unique_cell_type
names(marker_list)=names(COSG_markers[["names"]])


ref_com.FC_screen <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/ref_data/ref_com.FC_screen.rds")


for (i in 1:length(marker_list)) {
  marker_list[[i]]=marker_list[[i]][marker_list[[i]] %in% row.names(ref_com.FC_screen)]
}


saveRDS(marker_list,file='~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/inteRD/marker_list.rds')




marker_list <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/inteRD/marker_list.rds")
SCDC_output <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/scdc/SCDC_output.rds")


ref_data=SCDC_output[["prop.list"]]
for (i in 1:length(ref_data)) {
  ref_data[[i]]=as.matrix(ref_data[[i]][["prop.est.mvw"]])
}

unique_cell=c()
for (j in 1:length(ref_data)) {
  unique_cell=c(unique_cell,colnames(ref_data[[j]]))
}
unique_cell=unique(unique_cell)


for (k in 1:length(ref_data)) {
  add_celltype=unique_cell[!(unique_cell %in% colnames(ref_data[[k]]))]
  
  df_in=as.data.frame(ref_data[[k]])
  df_in[,add_celltype]<-0
  df_in=df_in[,sort(colnames(df_in))]
  ref_data[[k]]=as.matrix(df_in)
}


###构建marker_list和unique_cell_type
marker_list=marker_list[unique_cell]
unique_cell_type=unique_cell


# all(names(prop_results_list) == names(exp_scexp))
source("~/SWORD/interRD/interRD.R")



EXP <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/EXP_and_KEY/EXP_PBMC.rds")
######构建bulk
bulk.mtx<-as.data.frame(EXP[[1]])
for (i in 2:length(EXP)) {
  bulk.mtx2<-as.data.frame(EXP[[i]])
  bulk.mtx<-cbind(bulk.mtx,bulk.mtx2)
  
}
bulk.mtx<-as.matrix(bulk.mtx)
colnames(bulk.mtx)<-names(EXP)
fdata_bulk <- rownames(bulk.mtx)
pdata_bulk <- cbind(cellname = colnames(bulk.mtx), subjects = colnames(bulk.mtx))
bulk.mtx <- getESET(bulk.mtx, fdata = fdata_bulk, pdata = pdata_bulk)



lambda_option<-c(0,0.01,0.05,0.1,1,5,100)


# 
# ref_data[["Lee_SMC"]]<-NULL
# ref_data[["Khaliq"]]<-NULL

# unique_cell_type=colnames(ref_data[[1]])
# marker_list=marker_list[unique_cell]
library(Biobase)
source("~/SWORD/interRD/interRD.R")

InteRD1.output_in<-InteRD1lmh(bulk.data=bulk.mtx,
                              list_marker=marker_list,
                              cell_type_unique=unique_cell_type,
                              comb_used=ref_data,
                              lambda_option)

# save(bulk.mtx,marker_list,unique_cell_type,ref_data,file='~/ReCIDE/benchmark测试_最终/mult_reference_ER新/interRD/data2.rdata')


saveRDS(InteRD1.output_in,file = '~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/inteRD/interRD1_output.rds')

InteRD1_results<-InteRD.predict.prop(InteRD.output=InteRD1.output_in)
InteRD1_results=as.data.frame(t(InteRD1_results))

# key_df=readRDS("~/SWORD/单细胞图谱/参考集/整理/key_BRCA_ER_minor.rds")
# InteRD1_results=InteRD1_results[sort(row.names(InteRD1_results)),sort(colnames(InteRD1_results))]
# key_df=key_df[sort(row.names(key_df)),colnames(InteRD1_results)]
# 
# 
# sum(abs(InteRD1_results-key_df))/ncol(key_df)

saveRDS(InteRD1_results,file = '~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/inteRD/interRD1_prop_df.rds')
