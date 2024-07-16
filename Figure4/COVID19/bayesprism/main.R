t1= Sys.time ()

library(BayesPrism)
library(dplyr)

# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
# load.lbzip2('~/确定数据量影响/gse1760_26/ref_seurat.RDataFS',n.cores = 20)
library(fastSave) %>% suppressMessages() 
library(Seurat) %>% suppressMessages() 
library(pbmcapply) %>% suppressMessages() 
library(COSG) %>% suppressMessages() 


# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
# load.lbzip2('~/确定数据量影响/gse1760_26/ref_seurat.RDataFS',n.cores = 20)
ref_seurat=readRDS.lbzip2('~/scRNA_Seq_data/PBMC/Stephenson_seurat.rdsFS',n.cores = 200)
kk<-as.data.frame(ref_seurat@assays[["RNA"]]@counts)
kksum<-as.data.frame(apply(kk, 2, sum))
all(kksum>0)
kksum<-subset(kksum,kksum[,1]>0,)
kksum<-as.data.frame(kksum)
ref_seurat<-ref_seurat[,row.names(kksum)]
rm(kksum)
gc()

bulk.mtx<-readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/counts_data_query.rds")
rn=intersect(row.names(ref_seurat),row.names(bulk.mtx))
ref_seurat=ref_seurat[rn,]


bulk.mtx=as.data.frame(t(bulk.mtx))


###参考集构建

# #######################
###参考集构建
####celltype
cell.type.labels<-ref_seurat@meta.data['true']
cell.type.labels<-cell.type.labels[,1]

cell.type.labels=as.character(cell.type.labels)
cell.state.labels=NULL


######################

sc.dat.filtered.pc<-as.matrix(t(ref_seurat@assays[["RNA"]]@counts))

sc.dat.filtered.pc <- cleanup.genes (input=sc.dat.filtered.pc,
                                     input.type="count.matrix",
                                     species="hs", 
                                     gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                     exp.cells=5)
# plot.bulk.vs.sc (sc.input = sc.dat.filtered,
#                  bulk.input = bk.dat
#                  #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
# )

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

saveRDS(bp.res,file = '~/ReCIDE/应用_前二_新_inter/COVID19/bayesprism_benchmark/Bayes_output.rds')

saveRDS(theta,file = '~/ReCIDE/应用_前二_新_inter/COVID19/bayesprism_benchmark/Bayes_prd_df.rds')



Bayes_prd_df <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bayesprism_benchmark/Bayes_prd_df.rds")
# row.names(Bayes_prd_df)=names(EXP)
Bayes_prd_df=as.data.frame(t(Bayes_prd_df))
saveRDS(Bayes_prd_df,file="~/ReCIDE/应用_前二_新_inter/COVID19/bayesprism_benchmark/Bayes_prd_df.rds")


t2= Sys.time ()


