###参考的marker基因都是50个
###原始整理结果
###PVL直接原始结果
###CXCL13是细胞类型数大于40
sep_solDWLS <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/GSE164458/DWLS_output_ref50.rds")

for (i in 1:length(sep_solDWLS)) {
  
  sep_solDWLS[[i]]=sep_solDWLS[[i]][[1]]
  for (k in length(sep_solDWLS[[i]]):1) {
    # if(length(sep_solDWLS[[i]][[k]])<30){sep_solDWLS[[i]][[k]]=NULL}
    # if(!("PVL_Immature_s1" %in% names(sep_solDWLS[[i]][[k]]))){sep_solDWLS[[i]][[k]]=NULL}
    # if(!("PVL_Differentiated_s3" %in% names(sep_solDWLS[[i]][[k]]))){sep_solDWLS[[i]][[k]]=NULL}
      # if(!("T_cells_c3_CD4__Tfh_CXCL13" %in% names(sep_solDWLS[[i]][[k]]))){sep_solDWLS[[i]][[k]]=NULL}

  }
  
}

source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")

prd_after=ReCIDE_PCA(sep_solDWLS)

prd_after=prd_after[sort(row.names(prd_after)),sort(colnames(prd_after))]

saveRDS(prd_after,file='~/ReCIDE/应用_TNBC/论文中结果/GSE164458/prd_sep_原始.rds')



###参考的marker基因都是50个
###CXCL13_use
sep_solDWLS <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/GSE164458/DWLS_output_ref50.rds")

for (i in 1:length(sep_solDWLS)) {
  
  sep_solDWLS[[i]]=sep_solDWLS[[i]][[1]]
  for (k in length(sep_solDWLS[[i]]):1) {
    if(length(sep_solDWLS[[i]][[k]])<40){sep_solDWLS[[i]][[k]]=NULL}
    # if(!("PVL_Immature_s1" %in% names(sep_solDWLS[[i]][[k]]))){sep_solDWLS[[i]][[k]]=NULL}
    # if(!("PVL_Differentiated_s3" %in% names(sep_solDWLS[[i]][[k]]))){sep_solDWLS[[i]][[k]]=NULL}
    # if(!("T_cells_c3_CD4__Tfh_CXCL13" %in% names(sep_solDWLS[[i]][[k]]))){sep_solDWLS[[i]][[k]]=NULL}
    
  }
  
}

source("~/ReCIDE/PCA_and_hclust/ReCIDE_PCA.R")

prd_after=ReCIDE_PCA(sep_solDWLS)

prd_after=prd_after[sort(row.names(prd_after)),sort(colnames(prd_after))]

saveRDS(prd_after,file='~/ReCIDE/应用_TNBC/论文中结果/GSE164458/prd_sep_40.rds')

