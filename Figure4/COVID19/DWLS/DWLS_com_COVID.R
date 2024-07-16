library(DWLS)
library(pbmcapply)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(Seurat)


scdata_test<-readRDS.lbzip2('~/scRNA_Seq_data/PBMC/Stephenson_seurat.rdsFS',n.cores = 200)
bulkdata=readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/counts_data_query.rds")
rn=intersect(row.names(scdata_test),row.names(bulkdata))
scdata_test=scdata_test[rn,]

scdata_test@meta.data[["usetype"]]<-scdata_test@meta.data[["true"]]###########
cell.table<-as.data.frame(table(scdata_test@meta.data[["usetype"]]))
cell.table<-cell.table[cell.table[,2]>3,]
scdata_test<-subset(scdata_test,usetype %in% cell.table[,1])
gc()
labels<-as.character(scdata_test@meta.data[["usetype"]])
# labels <- str_replace(labels, "&", "_")
# labels <- str_replace(labels, "-", "_")
labels <- str_replace(labels, "-", "_")
labels <- str_replace(labels, "\\.", "_")
# labels <- str_replace(labels, " ", "_")

scdata_test@meta.data[['usetype']]<-labels
Idents(scdata_test)<-'usetype'

scdata_test<-NormalizeData(scdata_test)
# scdata_test<-NormalizeData(scdata_test)

COSG_markers <- cosg(
  scdata_test,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)

##整理marker基因名字等
de_group_list<-list()
for (i in 1:ncol(COSG_markers[["names"]])) {
  de_group_list[[i]]<-as.data.frame(COSG_markers[["names"]][,i])
  de_group_list[[i]][,2]<-as.data.frame(COSG_markers[["scores"]][,i])
  colnames(de_group_list[[i]])<-c('gene_name','scores')
  row.names(de_group_list[[i]])<-de_group_list[[i]][,1]
}
names(de_group_list)<-colnames(COSG_markers[["names"]])

##fold-change筛选
df2<-as.data.frame(t(as.data.frame(scdata_test@assays[["RNA"]]@data)))
gc()
df2[,"usetype"]<-scdata_test@meta.data[["usetype"]]
ref_list1<-group_by(df2,usetype) %>% summarize_each(mean)
ref_list1<-as.data.frame(ref_list1)
row.names(ref_list1)<-ref_list1[,1]
ref_list1<-ref_list1[,-1]
ref_list1<-as.data.frame(t(ref_list1))

secondFC <- c()
for(gene in rownames(ref_list1)){
  secondFC <- c(secondFC, sort(as.numeric(ref_list1[gene, ]), decreasing = T)[1]/sort(as.numeric(ref_list1[gene, ]), decreasing = T)[2])
}
names(secondFC) <- row.names(ref_list1)
secondFC<-secondFC[secondFC>1.5]
secondFC<-secondFC[!(is.na(secondFC))]

de_group_listFB<-de_group_list
for (i in 1:length(de_group_list)) {
  de_group_list[[i]]<-subset(de_group_list[[i]],row.names(de_group_list[[i]]) %in% names(secondFC))
}

#吃低保
# low_gene=10+round(length(de_group_listFB)/10)
# for (i in 1:length(de_group_list)) {
#   if(nrow(de_group_list[[i]])<10){
#     de_group_list[[i]]<-de_group_listFB[[i]][1:low_gene,]
#   }}



id<-labels
scdata<-as.matrix(scdata_test@assays[["RNA"]]@counts)
colnames(scdata) <- str_replace(colnames(scdata), "-", "_")


numberofGenes <- c()
for (i in unique(id)) {
  de_group <- de_group_list[[i]]
  # de_group <- readRDS(file = paste(path, "/de_", i, ".rds",
  #                                  sep = ""))
  # DEGenes <- rownames(de_group)[intersect(which(de_group$scores < pval.cutoff), which(de_group$avg_log2FC > diff.cutoff))]
  DEGenes <- rownames(de_group)
  
  nonMir = grep("MIR|Mir", DEGenes, invert = T)
  assign(paste("cluster_lrTest.table.", i, sep = ""), de_group[which(rownames(de_group) %in% DEGenes[nonMir]), ])
  numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))
}

f=100
Sig_list_in<-list()
conditionNumbers <- c()
for (G in 50:f) {
  Genes <- c()
  j = 1
  for (i in unique(id)) {
    if (numberofGenes[j] > 0) {
      temp <- paste("cluster_lrTest.table.", i, sep = "")
      temp <- as.name(temp)
      temp <- eval(parse(text = temp))
      temp <- temp[order(temp$scores, decreasing = TRUE),
      ]
      Genes <- c(Genes, (rownames(temp)[1:min(G, numberofGenes[j])]))
    }
    j = j + 1
  }
  Genes <- unique(Genes)
  ExprSubset <- scdata[Genes, ]
  Sig <- NULL
  # for (i in unique(id)) {
  #   Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id ==  i)]))))
  # }
  # 
  
  Sig<-as.data.frame(t(ExprSubset))
  Sig[,'id']<-id
  Sig=group_by(Sig,id) %>% summarize_each(mean)
  Sig<-as.data.frame(t(Sig))
  colnames(Sig)<-Sig[1,]
  Sig<-Sig[-1,]
  # Sig <- as.data.frame(lapply(Sig, as.numeric))
  Sig <- dplyr::mutate_all(Sig, as.numeric)
  
  # Sig2<-cbind(Sig,e1[row.names(Sig),])
  # kappa(Sig2)
  colnames(Sig) <- unique(id)
  conditionNumbers <- c(conditionNumbers, kappa(Sig))##降低kappa值，kappa值越高说明多重共线性越强，我们希望kappa值低一些
  Sig_list_in[[G-49]]<-Sig
}


# ##################
# plot(conditionNumbers)
# 
# ##每个细胞类型至少有10个marker，尝试让kappa值最小
# 
# # kappa<100则认为有共线性程度很小。
# # 100<=kappa<=1000则认为存在中等程度或者较强的共线性。
# # 若kappa>1000则认为存在很严重的共线性
G <-  which.min(conditionNumbers) + min(49, numberofGenes - 1)
#
#
Genes <- c()
j = 1
for (i in unique(id)) {
  if (numberofGenes[j] > 0) {
    temp <- paste("cluster_lrTest.table.", i, sep = "")
    temp <- as.name(temp)
    temp <- eval(parse(text = temp))
    temp <- temp[order(temp$scores, decreasing = TRUE),
    ]
    Genes <- c(Genes, (rownames(temp)[1:min(G, numberofGenes[j])]))
  }
  j = j + 1
}
Genes <- unique(Genes)
ExprSubset <- scdata[Genes, ]
# Sig <- NULL
# for (i in unique(id)) {
#   Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
# }
# colnames(Sig) <- unique(id)
Sig<-as.data.frame(t(ExprSubset))
Sig[,'id']<-id
Sig=group_by(Sig,id) %>% summarize_each(mean)
Sig<-as.data.frame(t(Sig))
colnames(Sig)<-Sig[1,]
Sig<-Sig[-1,]
# Sig <- as.data.frame(lapply(Sig, as.numeric))
Sig <- dplyr::mutate_all(Sig, as.numeric)


saveRDS(Sig,file = '~/ReCIDE/应用_前二_新_inter/COVID19/DWLS/PBMC_com_sig.rds')
# save(Sig_list,file = '~/SWORD/前聚类测试/CRC_Sig_list09.RData')



library(DWLS)
library(pbmcapply)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(Seurat)


EXP=readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/counts_data_query.rds")

Sig_list  <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/DWLS/PBMC_com_sig.rds")

###########################
com_ref_DWLS=Sig_list


sep_solDWLS<-list()
fun_DWLS_in<-function(i){
  bulk<-as.numeric(EXP[,i])
  # bulk=apply(bulk,2,sum)
  names(bulk)=row.names(EXP)
  
  
  ref1=as.matrix(com_ref_DWLS)
  
  tr<-trimData(ref1,bulk)
  output<-try(solveDampenedWLS(tr$sig,tr$bulk), TRUE)
  
  return(output)
}


sep_solDWLS<-pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores = 50)
names(sep_solDWLS)<-names(EXP)
# }
#save(sep_solDWLS)
# sep_solDWLS[[156]]<-NULL
saveRDS(sep_solDWLS,file='~/ReCIDE/应用_前二_新_inter/COVID19/DWLS/DWLS_com_output.rds')
# save(sep_solDWLS.list,file='~/SWORD/前聚类测试/sep_solDWLS_09ref.RData')

prd<-sep_solDWLS

prd[[1]]<-as.data.frame(prd[[1]])
df_merge<-prd[[1]]
colnames(df_merge)[1]<-names(prd)[1]

for(j in 2:length(prd)){
  
  prd[[j]]<-as.data.frame(prd[[j]])
  df_merge<-merge(df_merge, prd[[j]], by = "row.names", all = TRUE)
  row.names(df_merge)<-df_merge[,1]
  df_merge<-df_merge[,-1]
  colnames(df_merge)[j]<-names(prd)[j]
}

prd_com<-df_merge[,sort(names(df_merge))]

# prd_after['Mono_prolif',]=0
prd_com=prd_com[sort(row.names(prd_com)),sort(colnames(prd_com))]

saveRDS(prd_com,file='~/ReCIDE/应用_前二_新_inter/COVID19/DWLS/prd_com_PBMC.rds')


