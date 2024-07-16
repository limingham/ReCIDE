library(Seurat)
library(BayesPrism)
library(fastSave)
library(Seurat)
library(pbmcapply)
# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
combined.dgedata_150=readRDS.lbzip2('~/ReCIDE/benchmark测试_最终/high_res_TNBC/ref_data/ref_all.rdsFS',n.cores = 200)
combined.dgedata_150.list<- SplitObject(combined.dgedata_150, split.by = "batch")
combined.dgedata_150.list=combined.dgedata_150.list[sort(names(combined.dgedata_150.list))]

# ref_com.FC_screen <- readRDS("~/softlink/deone_pan/ref_data/ref_list.FC_screen.rds")

# all(names(combined.dgedata_150.list)==substr(names(ref_com.FC_screen),1,nchar(names(ref_com.FC_screen))-8))



func_makerPrismObject<-function(x){
  data_ref=subset(combined.dgedata_150,batch != names(combined.dgedata_150.list)[x])
  # data_ref=data_ref[row.names(ref_com.FC_screen[[x]]),]
  
  bulk_seurat<-subset(combined.dgedata_150,batch == names(combined.dgedata_150.list)[x])
  bulk.mtx<-as.matrix(apply(as.data.frame(bulk_seurat@assays[["RNA"]]@counts),1,sum))
  bulk.mtx<-cbind(bulk.mtx,bulk.mtx)
  bk.dat<-as.data.frame(t(bulk.mtx))
  # bk.dat=rbind(bk.dat,bk.dat)
  cell.type.labels<-data_ref@meta.data['celltype_subset']
  cell.type.labels<-cell.type.labels[,1]
  
  cell.type.labels<-as.character(cell.type.labels)
  cell.state.labels<-NULL
  
  sc.dat.filtered.pc<-as.matrix(t(data_ref@assays[["RNA"]]@counts))
  
  
  sc.dat.filtered.pc <- cleanup.genes (input=sc.dat.filtered.pc,
                                       input.type="count.matrix",
                                       species="hs", 
                                       gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                       exp.cells=5)
  
  
  
  myPrism <- new.prism(
    reference=sc.dat.filtered.pc, 
    mixture=bk.dat,
    input.type="count.matrix", 
    cell.type.labels = cell.type.labels, 
    cell.state.labels = cell.state.labels,
    key='Cancer Basal SC',
    outlier.cut=0.01,
    outlier.fraction=0.1)
  return(myPrism)
}
myPrism.list<-pbmclapply(1:length(combined.dgedata_150.list),func_makerPrismObject, mc.cores = 5)
names(myPrism.list)<-names(combined.dgedata_150.list)

saveRDS(myPrism.list,file='~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/TNBC_myPrism.list.rds')

myPrism.list=readRDS('~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/TNBC_myPrism.list.rds')



library(BayesPrism)
myPrism.list <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/TNBC_myPrism.list.rds")



theta<-list()
res.list<-list()
for (i in 1:length(myPrism.list)) {
  bp.res <- run.prism(prism = myPrism.list[[i]], n.cores=100)
  res.list[[i]]<-bp.res
  theta[[i]] <- get.fraction (bp=bp.res,
                              which.theta="final",
                              state.or.type="type")
  names(theta)[i]<-names(myPrism.list)[i]
}

BayesPrism.results<-theta
saveRDS(BayesPrism.results,file = '~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/BayesPrism_output.rds')


BayesPrism.results <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/BayesPrism_output.rds")


df_merge<-as.data.frame(t(BayesPrism.results[[1]])[,1])
colnames(df_merge)[1]<-names(BayesPrism.results)[1]

for(j in 2:length(BayesPrism.results)){
  
  BayesPrism.results[[j]]<-as.data.frame(t(BayesPrism.results[[j]])[,1])
  df_merge<-merge(df_merge, BayesPrism.results[[j]], by = "row.names", all = TRUE)
  row.names(df_merge)<-df_merge[,1]
  df_merge<-df_merge[,-1]
  colnames(df_merge)[j]<-names(BayesPrism.results)[j]
}
df_merge=df_merge[sort(row.names(df_merge)),sort(colnames(df_merge))]

prd_df=df_merge
saveRDS(prd_df,file = '~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/prd_df.rds')

# 
# 
# key_df <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/EXP_and_KEY/KEY.rds")
# prd_df <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/bayes_tutorial/prd_df.rds")
# prd_df[is.na(prd_df)]=0
# 
# prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
# key_df=key_df[sort(row.names(key_df)),sort(colnames(prd_df))]
# 
# RMSE_vec=c()
# for(i in 1:ncol(key_df)){
#   RMSE_vec[i]=ModelMetrics::rmse(as.numeric(key_df[,i]),as.numeric(prd_df[,i]))
#   # RMSE_vec[i]=cor(key_df[,i],prd_df[,i])
#   
# }
# mean(RMSE_vec)




