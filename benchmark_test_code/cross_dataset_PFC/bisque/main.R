library(Seurat)
library(BisqueRNA)
library(fastSave)
library(pbmcapply)


ref_seurat=readRDS.lbzip2(file='~/ReCIDE/benchmark测试/cross_dataset_PFC/ref_data/ref_all.rdsFS',n.cores = 200)
EXP<- readRDS("~/ReCIDE/benchmark测试/cross_dataset_PFC/EXP_and_key/EXP_PFC.rds")


##构建测试集
bulk.mtx<-as.data.frame(EXP[[1]])
for (i in 2:length(EXP)) {
  bulk.mtx2<-as.data.frame(EXP[[i]])
  bulk.mtx<-cbind(bulk.mtx,bulk.mtx2)
  
}
bulk.mtx<-as.matrix(bulk.mtx)
colnames(bulk.mtx)<-names(EXP)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.mtx)



###########################
kk<-as.data.frame(ref_seurat@assays[["RNA"]]@counts)
kksum<-as.data.frame(apply(kk, 2, sum))
all(kksum>0)
kksum<-subset(kksum,kksum[,1]>0,)
kksum<-as.data.frame(kksum)
ref_seurat<-ref_seurat[,row.names(kksum)]
rm(kksum)
rm(kk)

gc()

sample.ids <- colnames(ref_seurat)
individual.labels<-ref_seurat@meta.data['subject']
cell.type.labels<-ref_seurat@meta.data['major.celltype']

sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=individual.labels,
                       cellType=cell.type.labels)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)

sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(ref_seurat@assays$RNA@counts),
                                  phenoData=sc.pdata)
colnames(sc.eset@phenoData@data)<-c("SubjectName","cellType")
bisque_output <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)

# saveRDS(res,file = '~/SWORD/多参考集和单参考集比较测试/除SWORD外其他方法/bisque/TNBC_QEGAD_ref161529.rds')

saveRDS(bisque_output,file='~/ReCIDE/benchmark测试/cross_dataset_PFC/bisque/bisque_tutorial/bisque_output_new.rds')


prd_df <-as.data.frame(bisque_output[["bulk.props"]])
prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]


key_df<- readRDS("~/ReCIDE/benchmark测试/cross_dataset_PFC/EXP_and_key/key_PFC.rds")

library(dplyr)

prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]

# prd_df=prd_df[,c('P1','P2','P3','P4','P5','P6','P7','P8')]
# key_df=key_df[,c('P1','P2','P3','P4','P5','P6','P7','P8')]

rmse_vec=c()
for (i in 1:ncol(key_df)) {
  rmse_vec[i]=ModelMetrics::rmse(prd_df[,i],key_df[,i])
}

mean(rmse_vec)

saveRDS(prd_df,file='~/ReCIDE/benchmark测试/cross_dataset_PFC/bisque/bisque_tutorial/bisque_prd_df.rds')