library(DWLS)
library(pbmcapply)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(Seurat)


combine.dgedata_150=readRDS.lbzip2(file='~/ReCIDE/benchmark测试_最终/deone_kidney/ref_data/ref_all.rdsFS',n.cores = 200)
patient_names=unique(combine.dgedata_150@meta.data['patient'])


func_mu<-function(i){
  scdata_test<-subset(combine.dgedata_150,patient != patient_names[i,1])
  
  scdata_test@meta.data[["usetype"]]<-scdata_test@meta.data[["subclass.l1"]]###########
  cell.table<-as.data.frame(table(scdata_test@meta.data[["usetype"]]))
  cell.table<-cell.table[cell.table[,2]>3,]
  scdata_test<-subset(scdata_test,usetype %in% cell.table[,1])
  gc()
  labels<-as.character(scdata_test@meta.data[["usetype"]])
  labels <- str_replace(labels, "/", "_")
  # labels <- str_replace(labels, "-", "_")
  # labels <- str_replace(labels, " ", "_")
  # labels <- str_replace(labels, " ", "_")
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
  gene_names=c()
  for (i in 1:ncol(COSG_markers[["names"]])) {
    de_group_list[[i]]<-as.data.frame(COSG_markers[["names"]][,i])
    de_group_list[[i]][,2]<-as.data.frame(COSG_markers[["scores"]][,i])
    colnames(de_group_list[[i]])<-c('gene_name','scores')
    row.names(de_group_list[[i]])<-de_group_list[[i]][,1]
    gene_names=c(gene_names,COSG_markers[["names"]][[i]])
  }
  names(de_group_list)<-colnames(COSG_markers[["names"]])
  
  
  ############################
  ############################
  gene_names=unique(gene_names)
  scdata_test=scdata_test[gene_names,]
  ############################
  ############################
  
  ##fold-change筛选
  df2<-as.data.frame(t(as.data.frame(scdata_test@assays[["RNA"]]@data)))
  gc()
  df2[,"usetype"]<-scdata_test@meta.data[["usetype"]]
  ref_list1<-group_by(df2,usetype) %>% summarize_each(mean)
  ref_list1<-as.data.frame(ref_list1)
  row.names(ref_list1)<-ref_list1[,1]
  ref_list1<-ref_list1[,-1]
  ref_list1<-as.data.frame(t(ref_list1))
  
  rm(df2)
  # rm(scdata_test)
  
  gc()
  
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
  
  # 吃低保
  # low_gene=10+round(length(de_group_listFB)/10)
  # for (i in 1:length(de_group_list)) {
  #   if(nrow(de_group_list[[i]])<10){
  #     de_group_list[[i]]<-de_group_listFB[[i]][1:low_gene,]
  #   }}
  
  df3<-as.data.frame(t(as.data.frame(scdata_test@assays[["RNA"]]@counts)))
  gc()
  df3[,"usetype"]<-scdata_test@meta.data[["usetype"]]
  ref_mat<-group_by(df3,usetype) %>% summarize_each(mean)
  ref_mat<-as.data.frame(ref_mat)
  row.names(ref_mat)<-ref_mat[,1]
  ref_mat<-ref_mat[,-1]
  ref_mat<-as.data.frame(t(ref_mat))
  
  rm(df3)
  rm(scdata_test)
  gc()
  
  id<-labels
  # scdata<-as.matrix(scdata_test@assays[["RNA"]]@counts)
  # colnames(scdata) <- str_replace(colnames(scdata), "-", "_")
  
  
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
  
  f=200
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
    # ExprSubset <- scdata[Genes, ]
    # Sig <- NULL
    # for (i in unique(id)) {
    #   Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id ==  i)]))))
    # }
    # 
    
    # Sig<-as.data.frame(t(ExprSubset))
    # Sig[,'id']<-id
    # Sig=group_by(Sig,id) %>% summarize_each(mean)
    # Sig<-as.data.frame(t(Sig))
    # colnames(Sig)<-Sig[1,]
    # Sig<-Sig[-1,]
    # Sig <- dplyr::mutate_all(Sig, as.numeric)
    # 
    # Sig2<-cbind(Sig,e1[row.names(Sig),])
    # kappa(Sig2)
    # colnames(Sig) <- unique(id)
    Sig=ref_mat[Genes,]
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
  Sig=ref_mat[Genes,]
  Sig <- dplyr::mutate_all(Sig, as.numeric)
  
  
  if(nrow(Sig)<ncol(Sig)){
    conditionNumbers[which.min(conditionNumbers)] = conditionNumbers[which.min(conditionNumbers)]+1000
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
    Sig=ref_mat[Genes,]
    
  }
  
  
  Sig <- dplyr::mutate_all(Sig, as.numeric)
  
  return(Sig)
}
Sig_list=pbmclapply(1:nrow(patient_names),func_mu,mc.cores=5)
names(Sig_list)<-patient_names[,1]

saveRDS(Sig_list,file = '~/ReCIDE/benchmark测试_最终/deone_kidney/DWLS_com/kidney_com_sig.rds')
# save(Sig_list,file = '~/SWORD/前聚类测试/CRC_Sig_list09.RData')



# Sig_list<- readRDS("~/ReCIDE/benchmark测试/cross_dataset_PFC/ref_data/ref_com.FC_screen.rds")
Sig_list <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/DWLS_com/kidney_com_sig.rds")
Sig_list=Sig_list[sort(names(Sig_list))]


EXP<- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/EXP_and_KEY/EXP_kidney.rds")
EXP=EXP[names(Sig_list)]


###########################


sep_solDWLS<-list()
fun_DWLS_in<-function(i){
  bulk<-EXP[[i]]
  com_ref_DWLS=Sig_list[[i]]
  
  ref1=as.matrix(com_ref_DWLS)
  
  tr<-trimData(ref1,bulk)
  output<-try(solveDampenedWLS(tr$sig,tr$bulk), TRUE)
  
  return(output)
}


sep_solDWLS<-pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores = 10)
names(sep_solDWLS)<-names(EXP)
# }
#save(sep_solDWLS)
# sep_solDWLS[[156]]<-NULL
saveRDS(sep_solDWLS,file='~/ReCIDE/benchmark测试_最终/deone_kidney/DWLS_com/DWLS_com_output.rds')
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

# prd_com['NEU',]=0

saveRDS(prd_com,file='~/ReCIDE/benchmark测试_最终/deone_kidney/DWLS_com/prd_df.rds')




