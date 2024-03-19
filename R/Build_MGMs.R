
library(DWLS)
library(pbmcapply)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggbiplot)
library(pbmcapply)
library(mclust)
library(gmodels)


MGM_build<-function(SC_ref,EXP_df,label_celltype,label_subject,n_cores=5){
  n_cores_in=n_cores
  
  
  reference_1 <- ref_init(SC_ref,EXP_df,label_celltype,label_subject)
  ref_sep_FC_screen <- SecondFC_filter(reference_1,n_cores_in)
  
  
  ref_list<-Build_ReCIDE_Sig(reference_1,EXP_df,
                             ref_sep_FC_screen, n_cores_in)
  
  return(ref_list)
  
  
}


##############ref_build
##############
##############
ref_init<-function(SC_ref,EXP_df,label_celltype,label_subject){
  
  ####Left to reference and test shared genes and initially construct Seurat objects
  cell_names=intersect(row.names(SC_ref),row.names(EXP_df))
  SC_ref=SC_ref[cell_names,]
  
  if(class(SC_ref) != "Seurat"){
    SC_ref=CreateSeuratObject(SC_ref)
  }
  SC_ref@meta.data['celltype_label']=label_celltype
  SC_ref@meta.data['subject_label']=label_subject
  
  
  patient_id_df<-as.data.frame(table(SC_ref@meta.data[["subject_label"]]))
  patient_id_df<-patient_id_df[patient_id_df[,2]>500,]
  SC_ref<-subset(SC_ref,subject_label %in% patient_id_df[,1])
  
  cell_id=as.data.frame(table(SC_ref@meta.data[["celltype_label"]]))
  cell_id=cell_id[cell_id[,2]>20,]
  SC_ref<-subset(SC_ref,celltype_label %in% cell_id[,1])
  
  gc()
  
  
  
  ####100 marker genes were screened by cosg, a primary screen that ignores heterogeneity
  rowsum<-rowSums(SC_ref@assays[["RNA"]]@counts)
  rowsum<-as.data.frame(rowsum)
  rowsum2<-subset(rowsum,rowsum>0)
  SC_ref<-SC_ref[row.names(rowsum2),]
  SC_ref <- NormalizeData(SC_ref, verbose = FALSE)
  
  
  Idents(SC_ref)<-'celltype_label'
  COSG_markers <- cosg(
    SC_ref,
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
  
  
  
  # PANref.dgedata_150<-SC_ref[genelist[,1],]
  SC_ref.sub<-SC_ref[genelist[,1],]
  
  
  
  
  return(SC_ref.sub)
}


##############ref_build
##############
##############
SecondFC_filter <- function(SC_ref_sub,n_cores=5){
  
  patient_id<-as.data.frame(table(SC_ref_sub@meta.data[["subject_label"]]))
  patient_id<-patient_id[patient_id[,2]>500,]
  ##??????sign????????????
  func_in<-function(i){
    ref_seurat<-subset(SC_ref_sub,subject_label == as.character(patient_id[i,1]))
    ref_seurat<-NormalizeData(ref_seurat)
    df2<-as.data.frame(t(as.data.frame(ref_seurat@assays[["RNA"]]@data)))
    gc()
    df2[,"celltype_major"]<-ref_seurat@meta.data[["celltype_label"]]
    ref_list1<-group_by(df2,celltype_major) %>% summarize_each(mean)
    ref_list1<-as.data.frame(ref_list1)
    row.names(ref_list1)<-ref_list1[,1]
    ref_list1<-ref_list1[,-1]
    ref_list1<-as.data.frame(t(ref_list1))
    return(ref_list1)
  }
  
  ref_list<-pbmclapply(1:nrow(patient_id),func_in,mc.cores = n_cores)
  
  names(ref_list)=paste0(patient_id[,1],'_sep.ref')
  
  
  ref_listFB<-ref_list
  #ref_list<-ref_listFB
  
  
  
  func_secondFC<-function(i){
    secondFC <- c()
    for(gene in rownames(ref_list[[i]])){
      secondFC <- c(secondFC, sort(as.numeric(ref_list[[i]][gene, ]), decreasing = T)[1]/sort(as.numeric(ref_list[[i]][gene, ]), decreasing = T)[2])
    }
    names(secondFC) <- row.names(ref_list[[i]])
    # secondFC.list[[i]]<-secondFC
    # names(secondFC.list)[i]<-names(ref_list)[i]
    secondFC=sort(secondFC)
    return(secondFC)
  }
  
  secondFC.list<-pbmclapply(1:length(ref_list),func_secondFC,mc.cores = n_cores)
  names(secondFC.list)=paste0(patient_id[,1],'_sep.ref')
  
  func_findnames=function(i){
    ref_in=ref_list[[i]]
    n_col=ncol(ref_in)
    for (j in 1:nrow(ref_in)) {
      ref_in[j,n_col+1] = names(which.max(ref_in[j,]))
    }
    return(ref_in)
  }
  
  ref_list.names<-pbmclapply(1:length(ref_list),func_findnames,mc.cores = n_cores)
  names(ref_list.names)=paste0(patient_id[,1],'_sep.ref')
  
  
  
  func_FC_screen<-function(i){
    secondFC_in=secondFC.list[[i]][secondFC.list[[i]]>=1.5]
    ref_list.names_in=ref_list.names[[i]]
    colnames(ref_list.names_in)[ncol(ref_list.names_in)]='celltype'
    
    ref_list.names_in_FB=ref_list.names_in
    
    ref_list.names_in<-ref_list.names_in[names(secondFC_in),]
    
    gene_table=as.data.frame(table(ref_list.names_in[,ncol(ref_list.names_in)]))
    gene_table=gene_table[gene_table[,2]<10,]
    
    
    
    if(nrow(gene_table)>0){
      ref_list.names_in=subset(ref_list.names_in,!(celltype %in% gene_table[,1]))
      
      for(g in 1:nrow(gene_table)){
        cell_sub=ref_list.names_in_FB[ref_list.names_in_FB[,'celltype']==gene_table[g,1],]
        secondFC_sub=secondFC.list[[i]][row.names(cell_sub)]
        secondFC_sub=sort(secondFC_sub,decreasing = TRUE)
        if(length(secondFC_sub)>10){
          cell_sub=cell_sub[names(secondFC_sub)[1:10],]
        }else{cell_sub=cell_sub[names(secondFC_sub),]
        
        }
        ref_list.names_in=rbind(ref_list.names_in,cell_sub)
      }}
    
    ref_list.names_in=ref_list.names_in[,1:(ncol(ref_list.names_in)-1)]
    return(ref_list.names_in)
  }
  
  ref_list.FC_screen<-pbmclapply(1:length(ref_list.names),func_FC_screen,mc.cores = n_cores)
  
  names(ref_list.FC_screen)=paste0(patient_id[,1],'_sep.ref')
  
  return(ref_list.FC_screen)
}




##############ref_build
##############
##############

Build_ReCIDE_Sig<-function(SC_ref_sub,EXP_df,
                           ref_sep_names.FC_screen, n_cores=5){
  
  
  patient_id<-as.data.frame(table(SC_ref_sub@meta.data[["subject_label"]]))
  patient_id<-patient_id[patient_id[,2]>500,]
  
  
  ref_sep_names.FC_screen=ref_sep_names.FC_screen[sort(names(ref_sep_names.FC_screen))]
  
  
  rn=intersect(row.names(SC_ref_sub),names(EXP_df))
  SC_ref_sub=SC_ref_sub[rn,]
  
  
  
  sc_ref_list<-SplitObject(SC_ref_sub,split.by = 'subject_label')##############
  sc_ref_list=sc_ref_list[sort(names(sc_ref_list))]
  
  func_findf<-function(il){
    scdata_test<-sc_ref_list[[il]]
    scdata_test=scdata_test[row.names(ref_sep_names.FC_screen[[il]]),]
    
    scdata_test@meta.data[["usetype"]]<-scdata_test@meta.data[["celltype_label"]]###########
    cell.table<-as.data.frame(table(scdata_test@meta.data[["usetype"]]))
    cell.table<-cell.table[cell.table[,2]>3,]
    scdata_test<-subset(scdata_test,usetype %in% cell.table[,1])
    gc()
    labels<-as.character(scdata_test@meta.data[["usetype"]])
    # labels <- str_replace(labels, "&", "_")
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
    
    ##??????marker???????????????
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
    
    ##fold-change??????
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
    
    # ?????????
    low_gene=10+round(length(de_group_listFB)/10)
    for (i in 1:length(de_group_list)) {
      if(nrow(de_group_list[[i]])<10){
        de_group_list[[i]]<-de_group_listFB[[i]][1:low_gene,]
      }}
    
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
      conditionNumbers <- c(conditionNumbers, kappa(as.matrix(Sig)))##??????kappa??????kappa???????????????????????????????????????????????????kappa????????????
      Sig_list_in[[G-49]]<-Sig
    }
    
    
    # ##################
    # plot(conditionNumbers)
    # 
    # ##???????????????????????????10???marker????????????kappa?????????
    # 
    # # kappa<100????????????????????????????????????
    # # 100<=kappa<=1000??????????????????????????????????????????????????????
    # # ???kappa>1000????????????????????????????????????
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
  
  Sig_list<-pbmclapply(names(sc_ref_list),func_findf,mc.cores = n_cores)
  names(Sig_list)<-names(sc_ref_list)
  
  
  return(Sig_list)
}
