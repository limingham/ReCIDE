library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(dplyr)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 





SecondFC_filter <- function(SC_ref_sub,EXP_df,label_celltype,label_label_subject,n_cores=5){
  
  patient_id<-as.data.frame(table(SC_ref_sub@meta.data[["label_subject"]]))
  patient_id<-patient_id[patient_id[,2]>500,]
  ##制作sign基因矩阵
  func_in<-function(i){
    ref_seurat<-subset(SC_ref_sub,label_subject == as.character(patient_id[i,1]))
    ref_seurat<-NormalizeData(ref_seurat)
    df2<-as.data.frame(t(as.data.frame(ref_seurat@assays[["RNA"]]@data)))
    gc()
    df2[,"celltype_major"]<-ref_seurat@meta.data[["label_celltype"]]
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