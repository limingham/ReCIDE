library(dplyr)
library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 

##SC_ref: Reference counts matrix, can be in Seurat or dataframe format,rows represent genes,columns represent cells
#EXP_df：Targeted counts matrix, dataframe format,rows represent genes,columns represent cells
#label_celltype：Cell type labels, corresponding to the rows of SC_ref
#label_subject：Subject labels, corresponding to the rows of SC_ref

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
  
  
  
  
  ref_seurat<-SC_ref.sub
  ref_seurat<-NormalizeData(ref_seurat)
  df2<-as.data.frame(t(as.data.frame(ref_seurat@assays[["RNA"]]@data)))
  gc()
  df2[,"celltype_major"]<-ref_seurat@meta.data[["celltype_label"]]
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
  
  SC_ref_sub=SC_ref[row.names(ref_list.FC_screen),]
  
  return(SC_ref_sub)
}

