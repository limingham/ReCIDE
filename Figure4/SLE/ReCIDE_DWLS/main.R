library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(dplyr)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 


###先执行统一marker基因筛选
combined.data<-readRDS.lbzip2('~/scRNA_Seq_data/PBMC/Stephenson_seurat.rdsFS',n.cores = 200)
EXP_combine<- readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bulkdata/GSE50772_exprs.rds")

cell_names=intersect(row.names(combined.data),row.names(EXP_combine))
combined.data=combined.data[cell_names,]


patient_id_df<-as.data.frame(table(combined.data@meta.data[["patient_id"]]))
patient_id_df<-patient_id_df[patient_id_df[,2]>500,]
combined.data<-subset(combined.data,patient_id %in% patient_id_df[,1])

cell_id=as.data.frame(table(combined.data@meta.data[["true"]]))
cell_id=cell_id[cell_id[,2]>20,]
combined.data<-subset(combined.data,true %in% cell_id[,1])


gc()


rowsum<-rowSums(combined.data@assays[["RNA"]]@counts)
rowsum<-as.data.frame(rowsum)
rowsum2<-subset(rowsum,rowsum>0)
combined.data<-combined.data[row.names(rowsum2),]
combined.data <- NormalizeData(combined.data, verbose = FALSE)


Idents(combined.data)<-'true'
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

###再执行分患者marker基因筛选

saveRDS.lbzip2(combined.data.sub,file = '~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/ref_150.rdsFS',n.cores = 100)




library(fastSave) %>% suppressMessages() 
library(Seurat)  %>% suppressMessages() 
library(pbmcapply)  %>% suppressMessages() 
library(dplyr)  %>% suppressMessages() 
library(COSG)  %>% suppressMessages() 

combined.data.sub=readRDS.lbzip2('~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/ref_150.rdsFS',n.cores = 200)
# 
# bulk.mtx<- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/counts_data_query.rds")
# rn=intersect(row.names(combined.data.sub),row.names(bulk.mtx))
# combined.data.sub=combined.data.sub[rn,]

patient_id_df<-as.data.frame(table(combined.data.sub@meta.data[["patient_id"]]))
patient_id_df<-patient_id_df[patient_id_df[,2]>500,]
##制作sign基因矩阵
func_in<-function(i){
  ref_seurat<-subset(combined.data.sub,patient_id == as.character(patient_id_df[i,1]))
  ref_seurat<-NormalizeData(ref_seurat)
  df2<-as.data.frame(t(as.data.frame(ref_seurat@assays[["RNA"]]@data)))
  gc()
  df2[,"celltype_major"]<-ref_seurat@meta.data[["true"]]
  ref_list1<-group_by(df2,celltype_major) %>% summarize_each(mean)
  ref_list1<-as.data.frame(ref_list1)
  row.names(ref_list1)<-ref_list1[,1]
  ref_list1<-ref_list1[,-1]
  ref_list1<-as.data.frame(t(ref_list1))
  return(ref_list1)
}

ref_list<-pbmclapply(1:nrow(patient_id_df),func_in,mc.cores = 2)

names(ref_list)=paste0(patient_id_df[,1],'_sep.ref')


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

secondFC.list<-pbmclapply(1:length(ref_list),func_secondFC,mc.cores = 2)
names(secondFC.list)=paste0(patient_id_df[,1],'_sep.ref')

func_findnames=function(i){
  ref_in=ref_list[[i]]
  n_col=ncol(ref_in)
  for (j in 1:nrow(ref_in)) {
    ref_in[j,n_col+1] = names(which.max(ref_in[j,]))
  }
  return(ref_in)
}

ref_list.names<-pbmclapply(1:length(ref_list),func_findnames,mc.cores = 2)
names(ref_list.names)=paste0(patient_id_df[,1],'_sep.ref')



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

ref_list.FC_screen<-pbmclapply(1:length(ref_list.names),func_FC_screen,mc.cores = 2)

names(ref_list.FC_screen)=paste0(patient_id_df[,1],'_sep.ref')


saveRDS(ref_list.FC_screen, file="~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/ref_sep_names.FC_screen.rds")




####中间findmarker部分见PC端



library(DWLS)
library(pbmcapply)
library(fastSave)
library(COSG)
library(stringr)
library(dplyr)
library(Seurat)

EXP<<- readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bulkdata/GSE50772_exprs.rds")

# EXP=2^EXP-1
Sig_list<- readRDS("~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/PBMC_sep_sig.rds")
source('~/ReCIDE/ReCIDE_main/cosine_screen.R')
# EXP=EXP[unique(seurat_query_healthy@meta.data[["subject"]])]

sep_solDWLS<-list()
fun_DWLS_in<-function(i){
  bulk<-as.numeric(EXP[,i])
  # bulk=apply(bulk,2,sum)
  names(bulk)=row.names(EXP)
  sep_ref.list<-Sig_list
  
  
  sep_solDWLS_inside<-list()
  cos_list_inside=list()
  for (j in 1:length(sep_ref.list)) {
    ##dwls??????dataframe?????????matrix
    ref1<-as.matrix(sep_ref.list[[j]])
    
    # batch_output<-cosine_screen_HighToLow(ref1,bulk)
    batch_output<-cosine_screen_LowToHigh(ref1,bulk)
    # batch_output<-cosine_screen_GA(ref1,bulk)
    # 
    tr <- batch_output[[2]]
    tr<-trimData(ref1,bulk)
    sep_solDWLS_inside[[j]]<-try(solveDampenedWLS(tr[[1]],tr[[2]]), TRUE)
    
    cos_list_inside[[j]]<-batch_output[[1]]
    
    names(sep_solDWLS_inside)[j]<-names(sep_ref.list)[j]
    names(cos_list_inside)[j]<-names(sep_ref.list)[j]
    
    DWLS_output=list(sep_solDWLS_inside,cos_list_inside)
    names(DWLS_output)=c('sep_solDWLS','cos_list')
  }
  return(DWLS_output)
}


sep_solDWLS<-pbmclapply(1:length(EXP),fun_DWLS_in,mc.cores = 20)
names(sep_solDWLS)<-names(EXP)
# }

saveRDS(sep_solDWLS,file='~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/结果/DWLS_output_LowToHigh_PCA.rds')

sep_solDWLS<- readRDS("~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/结果/DWLS_output_LowToHigh_PCA.rds")
for (i in 1:length(sep_solDWLS)) {
  
  sep_solDWLS[[i]]=sep_solDWLS[[i]][[1]]
  
  
}

source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")

prd_after=ReCIDE_PCA(sep_solDWLS)

prd_after=prd_after[sort(row.names(prd_after)),sort(colnames(prd_after))]

saveRDS(prd_after,file='~/ReCIDE/应用_前二_新_inter/sep_ref_inter_SLE2/结果/prd_sep_LowToHigh_PCA.rds')



