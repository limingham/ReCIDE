library(BayesPrism)
library(dplyr)

# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
# load.lbzip2('~/确定数据量影响/gse1760_26/ref_seurat.RDataFS',n.cores = 20)
library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(dplyr)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 

ref_seurat<-readRDS.lbzip2(file='~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/ref_data/ref_all.rdsFS',n.cores = 200)
# kk<-as.data.frame(ref_seurat@assays[["RNA"]]@counts)
# kksum<-as.data.frame(apply(kk, 2, sum))
# all(kksum>0)
# kksum<-subset(kksum,kksum[,1]>0,)
# kksum<-as.data.frame(kksum)
# ref_seurat<-ref_seurat[,row.names(kksum)]
# rm(kksum)
# rm(kk)
# gc()


EXP <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/EXP_and_key/EXP_PFC.rds")
bulk.mtx<-as.data.frame(EXP[[1]])
for (i in 2:length(EXP)) {
  bulk.mtx2<-as.data.frame(EXP[[i]])
  bulk.mtx<-cbind(bulk.mtx,bulk.mtx2)
  
}
bulk.mtx<-as.data.frame(bulk.mtx)
bulk.mtx=as.data.frame(t(bulk.mtx))


###参考集构建
####celltype
cell.type.labels<-ref_seurat@meta.data['major.celltype']
cell.type.labels<-cell.type.labels[,1]

cell.type.labels=as.character(cell.type.labels)
cell.state.labels=NULL
# cell.state.labels<-paste0(cell.type.labels,'_',ref_seurat@meta.data[["subject"]])
# cell_table=as.data.frame(table(cell.state.labels))
# cell_table=cell_table[cell_table[,2]<10,]
# for (i in 1:length(cell.state.labels)) {
#   if(cell.state.labels[i] %in% cell_table[,1]){
#     cell.state.labels[i]=cell.type.labels[i]
#   }
# }  

######################

sc.dat.filtered.pc<-as.matrix(t(ref_seurat@assays[["RNA"]]@counts))

sc.dat.filtered.pc <- cleanup.genes (input=sc.dat.filtered.pc,
                                     input.type="count.matrix",
                                     species="hs", 
                                     gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                     exp.cells=5)




myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bulk.mtx,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1, )

bp.res <- run.prism(prism = myPrism, n.cores=100)
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")

saveRDS(bp.res,file = '~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/bayesprism/Bayes_output.rds')

saveRDS(theta,file = '~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/bayesprism/Bayes_prd_df.rds')



Bayes_prd_df <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/bayesprism/Bayes_prd_df.rds")
row.names(Bayes_prd_df)=names(EXP)
Bayes_prd_df=as.data.frame(t(Bayes_prd_df))
saveRDS(Bayes_prd_df,file="~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/bayesprism/Bayes_prd_df.rds")

