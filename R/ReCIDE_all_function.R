
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
# Build_ReCIDE_Sig(reference_1,EXP_df,
#                  ref_sep_FC_screen, n_cores_in)
Build_ReCIDE_Sig<-function(SC_ref_sub,EXP_df,
                           ref_sep_names.FC_screen, n_cores=5){
  
  
  patient_id<-as.data.frame(table(SC_ref_sub@meta.data[["subject_label"]]))
  patient_id<-patient_id[patient_id[,2]>500,]
  
  
  ref_sep_names.FC_screen=ref_sep_names.FC_screen[sort(names(ref_sep_names.FC_screen))]
  
  
  rn=intersect(row.names(SC_ref_sub),names(EXP_df))
  SC_ref_sub=SC_ref_sub[rn,]
  
  
  
  sc_ref_list<-SplitObject(SC_ref_sub,split.by = 'subject_label')##############
  sc_ref_list=sc_ref_list[sort(names(sc_ref_list))]
  
  for (l_in in length(sc_ref_list):1) {
    cell.table<-as.data.frame(table(sc_ref_list[[l_in]]@meta.data[["celltype_label"]]))
    cell.table<-cell.table[cell.table[,2]>3,]
    if(nrow(cell.table) <2){sc_ref_list = sc_ref_list[-l_in]}
    
    # if(length(unique(sc_ref_list[[l_in]]@meta.data[["celltype_label"]])) < 2){sc_ref_list = sc_ref_list[-l_in]}
    # print(length(unique(sc_ref_list[[l_in]]@meta.data[["celltype_label"]])))
  }
  
  
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
      conditionNumbers <- c(conditionNumbers, base::kappa(as.matrix(Sig)))##??????kappa??????kappa???????????????????????????????????????????????????kappa????????????
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






ReCIDE_deconvolution <- function(Sig_list,EXP_df,Method = 'DWLS',n_cores=5,n_celltype=0){
  
  EXP=list()
  for (i in 1:ncol(EXP_df)) {
    EXP[[i]]=as.numeric(EXP_df[,i])
    names(EXP[[i]])=row.names(EXP_df)
  }
  names(EXP)=colnames(EXP_df)
  
  
  
  if(Method == 'DWLS'){
    if(length(Sig_list)<1.5*length(EXP)){
      
      fun_DWLS_in<-function(i){
        bulk<-EXP[[i]]
        sep_ref.list<-Sig_list
        
        sep_solDWLS_inside<-list()
        for (j in 1:length(sep_ref.list)) {
          ##dwls??????dataframe?????????matrix
          ref1<-as.matrix(sep_ref.list[[j]])
          
          # batch_output<-cosine_screen_HighToLow(ref1,bulk)
          batch_output<-cosine_screen_LowToHigh(ref1,bulk)
          # batch_output<-cosine_screen_GA(ref1,bulk)
          
          tr <- batch_output[[2]]
          # tr<-trimData(ref1,bulk)
          sep_solDWLS_inside[[j]]<-try(solveDampenedWLS(tr[[1]],tr[[2]]), TRUE)
          
          
          names(sep_solDWLS_inside)[j]<-names(sep_ref.list)[j]
          
          
        }
        return(sep_solDWLS_inside)
      }
      sep_solDWLS=pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores=n_cores)
      names(sep_solDWLS)=names(EXP)
      
    }else{
      sep_solDWLS=list()
      
      
      for (f in 1:length(EXP)) {
        bulk<-EXP[[f]]
        sep_ref.list<-Sig_list
        
        DWLS_output<-list()
        
        
        fun_DWLS_in2<-function(j){
          ref1<-as.matrix(sep_ref.list[[j]])
          
          batch_output<-cosine_screen_LowToHigh(ref1,bulk)
          
          tr <- batch_output[[2]]
          sep_solDWLS_j<-try(solveDampenedWLS(tr[[1]],tr[[2]]), TRUE)
          
          return(sep_solDWLS_j)
          
        }
        sep_solDWLS_j=mclapply(1:length(sep_ref.list),fun_DWLS_in2,mc.cores=n_cores)
        names(sep_solDWLS_j)=names(sep_ref.list)
        
        sep_solDWLS[[f]]=sep_solDWLS_j
      }
      names(sep_solDWLS)=names(EXP)
      
    }
    
    
    for (j in 1:length(sep_solDWLS)) {
      for(k in length(sep_solDWLS[[j]]):1){
        if(length(sep_solDWLS[[j]][[k]])<n_celltype){sep_solDWLS[[j]][[k]]<-NULL}
      }
    }
    
  }else if(Method == 'CIBERSORT'){
    fun_DWLS_in<-function(i){
      bulk<-EXP[[i]]
      sep_ref.list<-Sig_list
      print(i)
      sep_solDWLS_inside<-list()
      for (j in 1:length(sep_ref.list)) {
        ##dwls??????dataframe?????????matrix
        ref1<-as.matrix(sep_ref.list[[j]])
        
        # batch_output<-cosine_screen_HighToLow(ref1,bulk)
        batch_output<-cosine_screen_LowToHigh(ref1,bulk)
        # batch_output<-cosine_screen_GA(ref1,bulk)
        
        tr <- batch_output[[2]]
        query_df=as.data.frame(tr[[2]])
        query_df=cbind(query_df,query_df)
        ref_df=as.data.frame(tr[[1]])
        
        sep_solDWLS_inside[[j]]<-CIBERSORT(ref_df,query_df,perm = 100,QN=TRUE)
        
        
        names(sep_solDWLS_inside)[j]<-names(sep_ref.list)[j]
        
        
      }
      return(sep_solDWLS_inside)
    }
    sep_solDWLS=pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores=n_cores)
    names(sep_solDWLS)=names(EXP)
    
  }
  

  prd_df=ReCIDE_PCA(sep_solDWLS,method=Method)
  
  prd_df=prd_df[sort(row.names(prd_df)),sort(colnames(prd_df))]
  
  
  deconvolution_output=list(sep_solDWLS,prd_df)
  names(deconvolution_output)=c('results_list','results_final_df')
  return(deconvolution_output)
}





ReCIDE_PCA=function(results_deconvolution,method='DWLS'){
  
  for (i in length(results_deconvolution):1) {
    for (j in length(results_deconvolution[[i]]):1) {
      if(class(results_deconvolution[[i]][[j]])[1]=="try-error"){results_deconvolution[[i]][[j]]<-NULL}
    }
  }
  
  if(method=='FARDEEP'){
    for (i in 1:length(results_deconvolution)) {
      for (j in 1:length(results_deconvolution[[i]])) {
        results_deconvolution[[i]][[j]]<-results_deconvolution[[i]][[j]][["relative.beta"]][1,]
      }}
    prd<-results_deconvolution
  }
  ##DWLS????????????
  # # prd<-sep_solOLS
  if(method=='DWLS'){
    for (i in length(results_deconvolution):1) {
      for (j in length(results_deconvolution[[i]]):1) {
        if(class(results_deconvolution[[i]][[j]])=="try-error"){results_deconvolution[[i]][[j]]<-NULL}
      }
    }
    prd<-results_deconvolution
  }
  
  
  # cibersort????????????
  if(method=='CIBERSORT'){
    for (i in 1:length(results_deconvolution)) {
      for (j in 1:length(results_deconvolution[[i]])) {
        data_in=results_deconvolution[[i]][[j]]
        data_in<-as.data.frame(data_in[1,-((ncol(data_in)-2):ncol(data_in))])
        data_in=as.numeric(data_in)
        names(data_in)=colnames(results_deconvolution[[i]][[j]])[1:(ncol(results_deconvolution[[i]][[j]])-3)]
        results_deconvolution[[i]][[j]]=data_in
      }
    }
    prd<-results_deconvolution
  }
  
  #prd??????????????????df??????
  prd_df<-list()
  for(i in 1:length(prd)){
    
    prd[[i]][[1]]<-as.data.frame(prd[[i]][[1]])
    df_merge<-prd[[i]][[1]]
    colnames(df_merge)[1]<-names(prd[[i]])[1]
    
    for(j in 2:length(prd[[i]])){
      
      prd[[i]][[j]]<-as.data.frame(prd[[i]][[j]])
      df_merge<-merge(df_merge, prd[[i]][[j]], by = "row.names", all = TRUE)
      row.names(df_merge)<-df_merge[,1]
      df_merge<-df_merge[,-1]
      colnames(df_merge)[j]<-names(prd[[i]])[j]
    }
    
    prd_df[[i]]<-as.data.frame(lapply(df_merge,as.numeric))
    row.names(prd_df[[i]])<-row.names(df_merge)
    names(prd_df)[i]<-names(prd)[i]
  }
  
  
  
  #prd_df<-prd_dfFB
  prd_dfFB<-prd_df
  for (j in 1:length(prd_df)) {
    for(i in 1:length(prd_df[[j]][,1])){
      ###0.75?????????????????????
      if(table(is.na(prd_df[[j]][i,]))['FALSE']<0.75*ncol(prd_df[[j]])){
        prd_df[[j]][i,][is.na(prd_df[[j]][i,])]<-0}
    }
  }
  ####?????????PCA????????????????????????PCA??????
  ##PCA
  ####PCA?????? 
  
  #################
  # for (n in 1:length(prd_df)) 
  func_pca<- function(n){
    
    data_train <- t(prd_df[[n]])
    
    data_train[is.na(data_train)]=0
    
    for (i in ncol(data_train):1) {
      #??????????????????0???????????????????????????????????????
      if((nrow(as.data.frame(data_train[data_train[,i]<=0,]))>0.5*nrow(data_train))){
        data_train=data_train[,-i]}
    }
    data_train=as.data.frame(data_train)
    if(ncol(data_train)>=3){
      library(gmodels)
      #????????????????????????????????????????????????
      PCA.OUT=fast.prcomp(data_train, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL)
      
      imp=PCA.OUT[["sdev"]]/sum(PCA.OUT[["sdev"]])
      
      
      pc_sure=2
      
      pca_n=PCA.OUT$x[,1:pc_sure]
      
      pca_n[,1]=scale(pca_n[,1])
      pca_n[,2]=scale(pca_n[,2])
      
      
      model <- Mclust(pca_n)
      
      
      # & (model_class[1,2]/nrow(pca_n)>0.5)
      ###????????????PCA??????????????????
      if(model$G !=1){
        model_class=as.data.frame(table(model$classification))
        model_class=model_class[order(model_class[,'Freq'],decreasing = TRUE),]
        
        if(model_class[1,2]>5){
          if(model_class[1,2]==model_class[2,2]){
            ##?????????????????????????????????
            
            model_class=model_class[model_class[,2]==model_class[1,2],]
            
            use_ind=model$classification[model$classification %in% model_class[,1]]
            data_train=data_train[names(use_ind),]
            
            data_train_list=split(data_train,use_ind) 
            
            kappa_value=c()
            for (ka in 1:length(data_train_list)) {
              data_in=as.data.frame(t(data_train_list[[ka]]))
              kappa_value[ka]=base::kappa(as.matrix(data_in))
            }
            
            
            patient_names=row.names(data_train_list[[which.max(kappa_value)]])
          }else{
            ##?????????????????????????????????
            patient_names=names(model$classification[model$classification==model_class[1,1]])
          }
        }else{
          ##??????else??????????????????5
          patient_names=colnames(prd_df[[n]])
        }
        
      }else{
        ##?????????????????????1???????????????
        patient_names=colnames(prd_df[[n]])}
    }
    
    return_data=prd_df[[n]][,patient_names]
    return(return_data)
  }
  
  prd_df2<-pbmclapply(1:length(prd_df),func_pca,mc.cores=100)
  names(prd_df2)<-names(prd_df)
  
  
  # for(p in 1:length(prd_df2)){
  #   if(ncol(prd_df2[[p]])>0.5*ncol(prd_df[[p]])){prd_df[[p]]=prd_df2[[p]]}
  # }
  prd_df=prd_df2
  prd_df_pca=prd_df2
  
  ####################################
  ###?????????????????????
  prd_mean<-list()
  for(i in 1:length(prd_df)){
    prd_mean[[i]]<-as.data.frame(prd_df[[i]][,1])
    row.names(prd_mean[[i]])<-row.names(prd_df[[i]])
    for(j in 1:length(prd_df[[i]][,1])){
      a1<-sort(as.numeric(prd_df[[i]][j,]))
      prd_mean[[i]][j,1]<-mean(as.numeric(a1))
      #prd_mean[[i]][j,1]<-mean(as.numeric(a1[,ceiling(length(a1)/5):floor(length(a1)/5*4)]))
    }
    names(prd_mean)[i]<-names(prd_df)[i]
    # prd_mean[[i]]<-prd_mean[[i]]/sum(prd_mean[[i]])
  }
  
  ##????????????????????????
  prd_median<-list()
  for(i in 1:length(prd_df)){
    prd_median[[i]]<-as.data.frame(prd_df[[i]][,1])
    row.names(prd_median[[i]])<-row.names(prd_df[[i]])
    for(j in 1:length(prd_df[[i]][,1])){
      a1<-sort(as.numeric(prd_df[[i]][j,]))
      prd_median[[i]][j,1]<-median(as.numeric(a1))
      #prd_mean[[i]][j,1]<-mean(as.numeric(a1[,ceiling(length(a1)/5):floor(length(a1)/5*4)]))
    }
    names(prd_median)[i]<-names(prd_df)[i]
    # prd_mean[[i]]<-prd_mean[[i]]/sum(prd_mean[[i]])
  }
  
  
  qua_sd<-list()
  qua_mean<-list()
  qua_noNA.list<-list()
  
  func_bagging<-function(j){
    print(paste('????????????',j,'???'))
    #df???bagging????????????????????????
    mean1<-prd_mean[[j]]
    mean2<-prd_mean[[j]]
    ##??????????????????
    sd2<-prd_mean[[j]]
    df1<-prd_df[[j]]
    qua_noNA<-list()
    for(i in 1:length(df1[,1])){
      ##Step2
      # df1[is.na(df1)]=0
      # noNA=df1[i,]
      noNA<-df1[i,which(!(is.na(df1[i,])))]
      
      
      noNA_sort<-sort(as.numeric(noNA))
      
      noNA=rbind(noNA,noNA)
      noNA<-noNA[,intersect(colnames(prd_df_pca[[j]]),colnames(noNA))]
      
      noNA=as.data.frame(t(noNA))
      
      ##Step3
      ##??????????????????????????????
      mean2[i,1]<-mean(as.numeric(noNA[,1]))
      sd2[i,1]<-sd(as.numeric(noNA[,1]))
      qua_noNA[[i]]<-noNA
      names(qua_noNA)[i]<-row.names(df1)[i]
      ##???mean1?????????mean2????????????????????????
    }
    
    return(list(sd2,mean2,qua_noNA))
  }
  
  all(names(prd_mean)==names(prd_df_pca))
  
  results_debagging<-pbmclapply(1:length(prd_mean),func_bagging,mc.cores=40)
  names(results_debagging)<-names(prd_mean)
  
  for (j in 1:length(results_debagging)) {
    qua_sd[[j]]<-results_debagging[[j]][[1]]
    qua_mean[[j]]<-results_debagging[[j]][[2]]
    qua_noNA.list[[j]]<-results_debagging[[j]][[3]]
    names(qua_mean)[j]<-names(results_debagging)[j]
    names(qua_noNA.list)[j]<-names(results_debagging)[j]
    names(qua_sd)[j]<-names(results_debagging)[j]
  }
  
  
  ##Step6
  ##????????????
  qua_noNA.listFB<-qua_noNA.list
  qua_meanFB<-qua_mean
  #qua_noNA.list<-qua_noNA.listFB
  for(j in 1:length(qua_noNA.listFB)){
    qua_meanFB[[j]][,1]<-qua_mean[[j]][,1]/sum(qua_mean[[j]][,1])
  }
  
  qua_meanFB_df<-qua_meanFB[[1]]
  colnames(qua_meanFB_df)[1]<-names(qua_meanFB)[1]
  
  for(j in 2:length(qua_meanFB)){
    
    qua_meanFB[[j]]<-as.data.frame(qua_meanFB[[j]])
    qua_meanFB_df<-merge(qua_meanFB_df, qua_meanFB[[j]], by = "row.names", all = TRUE)
    row.names(qua_meanFB_df)<-qua_meanFB_df[,1]
    qua_meanFB_df<-qua_meanFB_df[,-1]
    colnames(qua_meanFB_df)[j]<-names(qua_meanFB)[j]
  }
  qua_meanFB_df2<-qua_meanFB_df
  qua_meanFB_df<-as.data.frame(lapply(qua_meanFB_df2,as.numeric))
  row.names(qua_meanFB_df)<-row.names(qua_meanFB_df2)
  names(qua_meanFB_df)<-names(qua_meanFB)
  
  qua_meanFB_dfFB<-qua_meanFB_df
  for (j in 1:length(qua_meanFB_df[1,])) {
    qua_meanFB_df[,j][is.na(qua_meanFB_df[,j])]<-0}
  
  
  prd_after<-qua_meanFB_df[,sort(names(qua_meanFB_df))]
  
  return(prd_after)
}






cosine_screen_LowToHigh<-function(S,B){
  B=as.data.frame(B)
  B=apply(B,1,sum)
  set.seed(123)
  # B=B[B>0]
  S=S[intersect(names(B),row.names(S)),]
  B=B[intersect(names(B),row.names(S))]
  # 
  ##????????????second_FC??????????????????????????????????????????marker
  second_FC=c()
  for (i in 1:nrow(S)) {
    second_FC=c(second_FC,(sort(S[i,],decreasing = TRUE)[1]+0.01)/(sort(S[i,],decreasing = TRUE)[2]+0.01))
  }
  
  
  S2=as.data.frame(S)
  S2[,'secondFC_name']=names(second_FC)
  S2[,'secondFC_value']=second_FC
  
  unique_name=unique(names(second_FC))
  
  marker_screen=c()
  cos_list=list()
  for (n in 1:length(unique_name)) {
    
    ##???????????????marker???????????????
    second_FC_df_in=S2[S2[,'secondFC_name'] %in% unique_name[n],]
    second_FC_df_in[,'bulk_exp']=B[row.names(second_FC_df_in)]
    second_FC_df_in[,'gene_name']=row.names(second_FC_df_in)
    
    # second_FC_df_in=as.data.frame(second_FC_df_in)
    ##???secondFC??????????????????
    second_FC_df_in=arrange(second_FC_df_in,dplyr::desc(secondFC_value))
    row.names(second_FC_df_in)=second_FC_df_in[,'gene_name']
    
    
    ##?????????marker?????????????????????exp???
    mark_exp=second_FC_df_in[,c(unique_name[n],'bulk_exp')]
    mark_exp=mark_exp[mark_exp[,2]>0,]
    
    a=0
    b=0
    # cor1=cor(mark_exp[,1],mark_exp[,2],method='spearman')
    cor1=lsa::cosine(mark_exp[,1],mark_exp[,2])
    for (v in nrow(mark_exp):1) {
      if(nrow(mark_exp)<10){break}
      mark_exp_in=mark_exp[-v,]
      # cor_in=cor(mark_exp_in[,1],mark_exp_in[,2],method='spearman')
      cor_in=lsa::cosine(mark_exp_in[,1],mark_exp_in[,2])
      
      if(cor_in*nrow(mark_exp_in)>cor1*nrow(mark_exp)){
        mark_exp=mark_exp[-v,]
        cor1=cor_in
        # a=a+1
      }#else{b=b+1}
    }
    cos_list[[n]]=cor1
    marker_screen=c(marker_screen,row.names(mark_exp))
  }
  names(cos_list)=unique_name
  
  S_in=S[marker_screen,]
  data=as.data.frame(S_in)
  # data[,'query']=B[marker_screen]
  B_in=B[row.names(S_in)]
  
  list_cos=cos_list
  list_batch=list(as.matrix(S_in),B_in)
  return(list(list_cos,list_batch))
}



compare_group_proportion <- function(results_list,dir_diff_celltype_output){
  df_category1 = as.data.frame(t(results_list[[1]]))
  df_category2 = as.data.frame(t(results_list[[2]]))
  
  
  col_inter = intersect(colnames(df_category1),colnames(df_category2))
  diff1 <- setdiff(colnames(df_category1), col_inter)
  diff2 <- setdiff(colnames(df_category2), col_inter)
  if(length(diff1) > 0){df_category2[,diff1] = 0}
  
  if(length(diff2) > 0){df_category1[,diff2] = 0}
  
  df_category1 = df_category1[,sort(colnames(df_category1))]
  df_category2 = df_category2[,sort(colnames(df_category2))]
  
  p_values <- numeric(length = ncol(df_category1))
  trends <- character(length = ncol(df_category2))
  diff <- numeric(length = ncol(df_category2))
  
  for (j in seq_along(colnames(df_category1))) {
    col_name <- colnames(df_category1)[j]
    
    df_category1_in = df_category1[, col_name]
    
    
    df_category2_in = df_category2[, col_name]
    
    
    pvalue <- wilcox.test(df_category1_in, df_category2_in, alternative = "two.sided")$p.value
    p_values[j] <- pvalue
    diff[j] <- abs(median(df_category1_in)-median(df_category2_in))
    if (diff[j]==0){
      diff[j] <- abs(mean(df_category1_in)-mean(df_category2_in))
    }
    #   
    if (median(df_category1_in)<median(df_category2_in)) {
      trends[j] <- '+'  
    } else if (median(df_category1_in)>median(df_category2_in)) {
      trends[j] <- '-'
    } else {
      trends[j] <- '0'  
    }
  }
  adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
  
  re_df <- data.frame(
    Celltype = colnames(df_category1),
    PValue = adjusted_p_values,
    Trend = trends,
    Diff=diff
  )
  saveRDS(re_df,file=dir_diff_celltype_output)
}
