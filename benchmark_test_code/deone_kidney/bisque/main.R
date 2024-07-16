library(Seurat)
library(BisqueRNA)
library(fastSave)
library(pbmcapply)

# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
combine.dgedata_150=readRDS.lbzip2('~/ReCIDE/benchmark测试/deone_kidney/ref_data/ref_all.rdsFS',n.cores = 200)
patient_names=unique(combine.dgedata_150@meta.data['patient'])

bulk.mtx<-0
for (j in 1:nrow(patient_names)) {
  bulk_seurat<-subset(combine.dgedata_150,patient == patient_names[j,1])
  bulk.mtx1<-as.matrix(apply(as.data.frame(bulk_seurat@assays[["RNA"]]@counts),1,sum))
  if(length(bulk.mtx) != 1){
    bulk.mtx<-cbind(bulk.mtx,bulk.mtx1)}else{
      bulk.mtx<-bulk.mtx1}
}

bulk.mtx<-as.matrix(bulk.mtx)
colnames(bulk.mtx)<-patient_names[,1]
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.mtx)


##############################################
function_bisque_in<-function(i){
  sc_seurat<-subset(combine.dgedata_150,patient != patient_names[i,1])
  # sc_seurat=sc_seurat[row.names(ref_com.FC_screen[[i]]),]
  cell_num=apply(sc_seurat@assays$RNA@counts,2,sum)
  cell_num=cell_num[cell_num>0]
  sc_seurat=sc_seurat[,names(cell_num)]
  gc()
  sample.ids <- colnames(sc_seurat)
  individual.labels<-sc_seurat@meta.data['patient']
  cell.type.labels<-sc_seurat@meta.data['subclass.l1']
  
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
  
  sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(sc_seurat@assays$RNA@counts),
                                    phenoData=sc.pdata)
  gc()
  
  colnames(sc.eset@phenoData@data)<-c("SubjectName","cellType")
  
  
  res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)
  results<-as.data.frame(res[["bulk.props"]][,as.character(patient_names[i,1])])
  colnames(results)<-patient_names[i,1]
  return(results)
}
results_bisque<-pbmclapply(1:nrow(patient_names),function_bisque_in,mc.cores = 3)

names(results_bisque)<-patient_names[,1]

saveRDS(results_bisque,file = '~/ReCIDE/benchmark测试/deone_kidney/bisque/bisque_tutorial/bisuqe_output.rds')

prd<-results_bisque

df_merge<-as.data.frame(prd[[1]])
# colnames(df_merge)[1]<-names(prd)[1]
colnames(df_merge)[1]<-names(prd[[1]])

for(j in 2:length(prd)){
  
  prd[[j]]<-as.data.frame(prd[[j]])
  df_merge<-merge(df_merge, prd[[j]], by = "row.names", all = TRUE)
  row.names(df_merge)<-df_merge[,1]
  df_merge<-df_merge[,-1]
  # colnames(df_merge)[j]<-names(prd)[j]
  colnames(df_merge)[j]<-names(prd[[j]])
}



df_merge<-as.data.frame(df_merge)
df_merge<-df_merge[,sort(colnames(df_merge))]
df_merge<-df_merge[sort(row.names(df_merge)),]


saveRDS(df_merge,file='~/ReCIDE/benchmark测试/deone_kidney/bisque/bisque_tutorial/prd_bisque.rds')
prd_df=df_merge
key_df <- readRDS("~/ReCIDE/benchmark测试/deone_kidney/EXP_and_KEY/KEY_kidney.rds")

prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
key_df=key_df[sort(row.names(key_df)),sort(colnames(prd_df))]

RMSE_vec=c()
for(i in 1:ncol(key_df)){
  RMSE_vec[i]=ModelMetrics::rmse(as.numeric(key_df[,i]),as.numeric(prd_df[,i]))
  # RMSE_vec[i]=cor(key_df[,i],prd_df[,i])
  
}
median(RMSE_vec)
