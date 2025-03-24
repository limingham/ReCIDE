library(DWLS)
library(Seurat)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(pbmcapply)
library(mclust)
library(gmodels)

func_run_ReCIDE = function(SC_ref,EXP_df,celltype,subject,dir_ref,dir_results){
  
  celltype_label = SC_ref@meta.data[[celltype]]
  subject_label = SC_ref@meta.data[[subject]]
  
  ref_MGM.list = MGM_build(SC_ref = SC_ref,
                           EXP_df = EXP_df,
                           label_celltype = celltype_label,
                           label_subject = subject_label,
                           n_cores = 100)
  
  saveRDS(ref_MGM.list,file = dir_ref)
  
  
  # ref_MGM.list <- readRDS("~/ReCIDE/benchmark_syq/ReCIDE/normal_ref2.rds")
  
  ReCIDE_results = ReCIDE_deconvolution(Sig_list = ref_MGM.list,
                                        EXP_df = EXP_df,
                                        Method = 'DWLS',
                                        n_cores = 100,
                                        n_celltype = 0)
  
  
  
  saveRDS(ReCIDE_results,file = dir_results)
}



func_run_ReCIDE_and_postprocessing = function(SC_ref1,SC_ref2,EXP_df1,EXP_df2,celltype1,celltype2,subject1,subject2,dir_ref1,dir_ref2,dir_results1,dir_results2,dir_diff_celltype_output){
  func_run_ReCIDE(SC_ref1,EXP_df1,celltype1,subject1,dir_ref1,dir_results1)
  func_run_ReCIDE(SC_ref2,EXP_df2,celltype2,subject2,dir_ref2,dir_results2)
  
  results_list=list()
  results_list[[1]]=readRDS(dir_results1)[["results_final_df"]]
  results_list[[2]]=readRDS(dir_results2)[["results_final_df"]]
  
  compare_group_proportion(results_list,dir_diff_celltype_output)
}





library(Seurat)
library(BisqueRNA)
library(fastSave)
library(pbmcapply)

func_run_bisque = function(ref_seurat,bulkdata,celltype,subject,dir_results){
  bulk.mtx=as.matrix(bulkdata)
  
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.mtx)
  ###########################
  kk<-as.data.frame(ref_seurat@assays[["RNA"]]@counts)
  kksum<-as.data.frame(apply(kk, 2, sum))
  all(kksum>0)
  kksum<-subset(kksum,kksum[,1]>0,)
  kksum<-as.data.frame(kksum)
  ref_seurat<-ref_seurat[,row.names(kksum)]
  rm(kksum)
  gc()
  
  sample.ids <- colnames(ref_seurat)
  individual.labels<-ref_seurat@meta.data[[subject]]
  cell.type.labels<-ref_seurat@meta.data[[celltype]]
  
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
  
  
  prd_df <-as.data.frame(bisque_output[["bulk.props"]])
  prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
  saveRDS(prd_df,file=dir_results)
}



library(dplyr)
library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(dplyr)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 


func_run_CIBERSORT = function(combined.data,bulk.mtx,celltype,subject,dir_ref,dir_results){
  
  combined.data@meta.data['celltype_label']=combined.data@meta.data[[celltype]]
  combined.data@meta.data['subject_label']=combined.data@meta.data[[subject]]
  
  
  patient_id_df<-as.data.frame(table(combined.data@meta.data[["subject_label"]]))
  patient_id_df<-patient_id_df[patient_id_df[,2]>500,]
  combined.data<-subset(combined.data,subject_label %in% patient_id_df[,1])
  
  cell_id=as.data.frame(table(combined.data@meta.data[["celltype_label"]]))
  cell_id=cell_id[cell_id[,2]>20,]
  combined.data<-subset(combined.data,celltype_label %in% cell_id[,1])
  
  gc()
  
  
  rowsum<-rowSums(combined.data@assays[["RNA"]]@counts)
  rowsum<-as.data.frame(rowsum)
  rowsum2<-subset(rowsum,rowsum>0)
  combined.data<-combined.data[row.names(rowsum2),]
  combined.data <- NormalizeData(combined.data, verbose = FALSE)
  
  
  Idents(combined.data)<- celltype
  COSG_markers <- cosg(
    combined.data,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  
  
  genelist<-as.data.frame(COSG_markers[["names"]][,1])
  colnames(genelist)<-'dge'
  for(gene in 2:ncol(COSG_markers[["names"]])){
    dge2<-as.data.frame(COSG_markers[["names"]][,gene])
    colnames(dge2)<-'dge'
    genelist<-rbind(genelist,dge2)
  }
  genelist<-unique(genelist)
  
  
  
  # PANref.dgedata_150<-combined.data[genelist[,1],]
  combined.data.sub<-combined.data[genelist[,1],]
  
  #saveRDS.lbzip2(combined.data.sub, file="~/ReCIDE/benchmark??????_??????/cross_dataset_PFC/ref_data/PFC_cross.dgedata_150.rdsFS",n.cores = 200)
  
  
  
  ref_seurat<-combined.data.sub
  ref_seurat<-NormalizeData(ref_seurat)
  df2<-as.data.frame(t(as.data.frame(ref_seurat@assays[["RNA"]]@data)))
  gc()
  df2[,"celltype_major"]<-ref_seurat@meta.data[[celltype]]
  ref_list<-group_by(df2,celltype_major) %>% summarize_each(mean)
  ref_list<-as.data.frame(ref_list)
  row.names(ref_list)<-ref_list[,1]
  ref_list<-ref_list[,-1]
  ref_list<-as.data.frame(t(ref_list))
  
  
  
  ref_listFB<-ref_list
  #ref_list<-ref_listFB
  
  secondFC <- c()
  for(gene in rownames(ref_list)){
    secondFC <- c(secondFC, sort(as.numeric(ref_list[gene, ]), decreasing = T)[1]/sort(as.numeric(ref_list[gene, ]), decreasing = T)[2])
  }
  names(secondFC) <- row.names(ref_list)
  # secondFC<-secondFC
  # names(secondFC.list)[i]<-names(ref_list)[i]
  secondFC=sort(secondFC)
  
  
  
  ref_in=ref_list
  n_col=ncol(ref_in)
  for (j in 1:nrow(ref_in)) {
    ref_in[j,n_col+1] = names(which.max(ref_in[j,]))
  }
  
  
  
  
  secondFC_in=secondFC[secondFC>=1.5]
  ref_list.names_in=ref_in
  colnames(ref_list.names_in)[ncol(ref_list.names_in)]='celltype'
  
  ref_list.names_in_FB=ref_list.names_in
  
  ref_list.names_in<-ref_list.names_in[names(secondFC_in),]
  
  gene_table=as.data.frame(table(ref_list.names_in[,ncol(ref_list.names_in)]))
  gene_table=gene_table[gene_table[,2]<10,]
  
  
  
  if(nrow(gene_table)>0){
    ref_list.names_in=subset(ref_list.names_in,!(celltype %in% gene_table[,1]))
    
    for(g in 1:nrow(gene_table)){
      cell_sub=ref_list.names_in_FB[ref_list.names_in_FB[,'celltype']==gene_table[g,1],]
      secondFC_sub=secondFC[row.names(cell_sub)]
      secondFC_sub=sort(secondFC_sub,decreasing = TRUE)
      if(length(secondFC_sub)>10){
        cell_sub=cell_sub[names(secondFC_sub)[1:10],]
      }else{cell_sub=cell_sub[names(secondFC_sub),]
      
      }
      ref_list.names_in=rbind(ref_list.names_in,cell_sub)
    }}
  
  ref_list.FC_screen=ref_list.names_in[,1:(ncol(ref_list.names_in)-1)]
  
  
  
  
  
  saveRDS(ref_list.FC_screen, file=dir_ref)
  
  
  CoreAlg <- function(X, y){
    
    #try different values of nu
    svn_itor <- 3
    
    res <- function(i){
      if(i==1){nus <- 0.25}
      if(i==2){nus <- 0.5}
      if(i==3){nus <- 0.75}
      model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
      model
    }
    
    if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
      out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)
    
    nusvm <- rep(0,svn_itor)
    corrv <- rep(0,svn_itor)
    
    #do cibersort
    t <- 1
    while(t <= svn_itor) {
      weights = t(out[[t]]$coefs) %*% out[[t]]$SV
      weights[which(weights<0)]<-0
      w<-weights/sum(weights)
      u <- sweep(X,MARGIN=2,w,'*')
      k <- apply(u, 1, sum)
      nusvm[t] <- sqrt((mean((k - y)^2)))
      corrv[t] <- cor(k, y)
      t <- t + 1
    }
    
    #pick best model
    rmses <- nusvm
    mn <- which.min(rmses)
    model <- out[[mn]]
    
    #get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q<0)]<-0
    w <- (q/sum(q))
    
    mix_rmse <- rmses[mn]
    mix_r <- corrv[mn]
    
    newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
    
  }
  
  #' do permutations
  #' @param perm Number of permutations
  #' @param X cell-specific gene expression
  #' @param y mixed expression per sample
  #' @export
  doPerm <- function(perm, X, Y){
    itor <- 1
    Ylist <- as.list(data.matrix(Y))
    dist <- matrix()
    
    while(itor <= perm){
      #print(itor)
      
      #random mixture
      yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
      
      #standardize mixture
      yr <- (yr - mean(yr)) / sd(yr)
      
      #run CIBERSORT core algorithm
      result <- CoreAlg(X, yr)
      
      mix_r <- result$mix_r
      
      #store correlation
      if(itor == 1) {dist <- mix_r}
      else {dist <- rbind(dist, mix_r)}
      
      itor <- itor + 1
      if(itor %% 50 == 0){
        print(paste('????????????doPerm??????',itor,'?????????'))}
    }
    newList <- list("dist" = dist)
  }
  
  #' Main functions
  #' @param sig_matrix file path to gene expression from isolated cells
  #' @param mixture_file heterogenous mixed expression
  #' @param perm Number of permutations
  #' @param QN Perform quantile normalization or not (TRUE/FALSE)
  #' @export
  CIBERSORT <- function(sig_matrix, mixture_file, perm=100, QN=TRUE){
    
    #read in data
    # X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
    # Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)
    
    X <- data.matrix(sig_matrix)
    Y <- data.matrix(mixture_file)
    
    #order
    X <- X[order(rownames(X)),]
    Y <- Y[order(rownames(Y)),]
    
    P <- perm #number of permutations
    
    #anti-log if max < 50 in mixture file
    if(max(Y) < 50) {Y <- 2^Y}
    
    #quantile normalization of mixture file
    if(QN == TRUE){
      tmpc <- colnames(Y)
      tmpr <- rownames(Y)
      Y <- preprocessCore::normalize.quantiles(Y)
      colnames(Y) <- tmpc
      rownames(Y) <- tmpr
    }
    
    #intersect genes
    Xgns <- row.names(X)
    Ygns <- row.names(Y)
    YintX <- Ygns %in% Xgns
    Y <- Y[YintX,]
    XintY <- Xgns %in% row.names(Y)
    X <- X[XintY,]
    
    #standardize sig matrix
    X <- (X - mean(X)) / sd(as.vector(X))
    
    #empirical null distribution of correlation coefficients
    if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
    
    print(nulldist)
    
    header <- c("Sample",colnames(X),"P-value","Correlation","RMSE")
    #header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
    #print(header)
    
    output <- matrix()
    itor <- 1
    mixtures <- dim(Y)[2]
    pval <- 9999
    
    #iterate through mixtures
    while(itor <= mixtures){
      
      y <- Y[,itor]
      
      #standardize mixture
      y <- (y - mean(y)) / sd(y)
      
      #run SVR core algorithm
      result <- CoreAlg(X, y)
      
      #get results
      w <- result$w
      mix_r <- result$mix_r
      mix_rmse <- result$mix_rmse
      
      #calculate p-value
      if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
      
      #print output
      out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
      if(itor == 1) {output <- out}
      else {output <- rbind(output, out)}
      
      itor <- itor + 1
      
    }
    
    #save results
    write.table(rbind(header,output), file= dir_results, sep="\t", row.names=F, col.names=F, quote=F)
    
    #return matrix object containing all results
    # row.names(output)=output[,1]
    # output<-as.data.frame(output[,-1])
    # obj <- rbind(header,output)
    # 
    # obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
    #colnames(output) <- c(colnames(X),"P-value","Correlation","RMSE")
    #return(output)
  }
  sig_matrix <- readRDS(dir_ref)
  mixture_file <- bulk.mtx
  CIBERSORT(sig_matrix,mixture_file)
}




suppressMessages(library(MuSiC))
suppressMessages(library(Seurat))
suppressMessages(library(fastSave))
suppressMessages(library(pbmcapply))
suppressMessages(library(SingleCellExperiment))

func_run_DWLS = function(scdata_test,bulkdata,celltype,dir_DWLS_ref,dir_DWLS_results){
  rn=intersect(row.names(scdata_test),row.names(bulkdata))
  scdata_test=scdata_test[rn,]
  
  
  scdata_test@meta.data[["usetype"]]<-scdata_test@meta.data[[celltype]]
  cell.table<-as.data.frame(table(scdata_test@meta.data[["usetype"]]))
  cell.table<-cell.table[cell.table[,2]>3,]
  scdata_test<-subset(scdata_test,usetype %in% cell.table[,1])
  gc()
  labels<-as.character(scdata_test@meta.data[["usetype"]])
  
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
    print(G)
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
  
  
  saveRDS(Sig,file = dir_DWLS_ref)
  
  
  # EXP<- readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bulkdata/GSE50772_exprs.rds")
  
  Sig_list <- readRDS(dir_DWLS_ref)
  
  ###########################
  com_ref_DWLS=Sig_list
  EXP = bulkdata
  
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
  
  saveRDS(prd_com,file = dir_DWLS_results)
}





suppressMessages(library(MuSiC))
suppressMessages(library(Seurat))
suppressMessages(library(fastSave))
suppressMessages(library(pbmcapply))
suppressMessages(library(SingleCellExperiment))


func_run_music = function(ref_seurat,EXP_df,celltype,subject,dir_results){
  
  bulk.mtx<-as.matrix(EXP_df)
  
  
  kk<-as.data.frame(ref_seurat@assays[["RNA"]]@counts)
  kksum<-as.data.frame(apply(kk, 2, sum))
  all(kksum>0)
  kksum<-subset(kksum,kksum[,1]>0,)
  kksum<-as.data.frame(kksum)
  ref_seurat<-ref_seurat[,row.names(kksum)]
  rm(kksum)
  
  
  gc()
  
  
  #单细胞参考集是sce的形式
  sce <- SingleCellExperiment(as.matrix(ref_seurat@assays$RNA@counts),
                              colData=DataFrame(label=ref_seurat@meta.data),
                              rowData=DataFrame(length=row.names(ref_seurat)))
  names(assays(sce))<-'counts'
  
  #在music前要先确定bulk和sc是否只有一列，只有一列的话无法运行（包括sample和cluster数是否为1）
  music_output = music_prop(bulk.mtx = bulk.mtx, sc.sce = sce, clusters = paste0('label.',celltype1),
                            samples = paste0('label.',subject1),verbose = F)
  
  
  saveRDS(music_output,file= dir_results)
}



library(BayesPrism)
library(dplyr)

# load("~/代表性方法测试/BayesPrism/tutorial/tutorial.gbm.rdata")
# load.lbzip2('~/确定数据量影响/gse1760_26/ref_seurat.RDataFS',n.cores = 20)
library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 




func_run_bayes = function(ref_seurat,bulk.mtx,celltype,dir_results,key1=NULL){
  kk<-as.data.frame(ref_seurat@assays[["RNA"]]@counts)
  kksum<-as.data.frame(apply(kk, 2, sum))
  all(kksum>0)
  kksum<-subset(kksum,kksum[,1]>0,)
  kksum<-as.data.frame(kksum)
  ref_seurat<-ref_seurat[,row.names(kksum)]
  rm(kksum)
  gc()
  
  # bulk.mtx = as.data.frame(readRDS("~/ReCIDE/benchmark_syq/CRC/bulk/counts_GSE39582_dMMR.rds"))
  rn=intersect(row.names(ref_seurat),row.names(bulk.mtx))
  ref_seurat=ref_seurat[rn,]
  
  
  bulk.mtx=as.data.frame(t(bulk.mtx))
  
  
  ###参考集构建
  
  # #######################
  ###参考集构建
  ####celltype
  cell.type.labels<-ref_seurat@meta.data[[celltype]]
  # cell.type.labels<-cell.type.labels[,1]
  
  cell.type.labels=as.character(cell.type.labels)
  cell.state.labels=NULL
  
  
  ######################
  
  sc.dat.filtered.pc<-as.matrix(t(ref_seurat@assays[["RNA"]]@counts))
  
  sc.dat.filtered.pc <- cleanup.genes(input=sc.dat.filtered.pc,
                                      input.type="count.matrix",
                                      species="hs", 
                                      gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                      exp.cells=5)
  # plot.bulk.vs.sc (sc.input = sc.dat.filtered,
  #                             bulk.input = bk.dat
  #                             #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
  # )
  
  myPrism <- new.prism(
    reference=sc.dat.filtered.pc, 
    mixture=bulk.mtx,
    input.type="count.matrix", 
    cell.type.labels = cell.type.labels, 
    cell.state.labels = cell.state.labels,
    key=key1,
    outlier.cut=0.01,
    outlier.fraction=0.1, )
  
  bp.res <- run.prism(prism = myPrism, n.cores=50)
  theta <- get.fraction (bp=bp.res,
                         which.theta="final",
                         state.or.type="type")
  
  
  
  saveRDS(theta,file = dir_results)
}
