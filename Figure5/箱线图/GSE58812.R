library(survival)
library(survminer)
prd_after <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/GSE58812/prd_sep_init.rds")
Clinical<- readRDS("~/ReCIDE/应用_TNBC/GSE58812/GSE58812_metadata.rds")

pvl_all=prd_after['PVL_Differentiated_s3',]+prd_after['PVL_Immature_s1',]+prd_after['PVL_Immature_s2',]+prd_after['Cycling_PVL',]

data_matrix<-data.frame(CXCL13 = as.numeric(prd_after['T_cells_c3_CD4__Tfh_CXCL13',]),
                        PVLs1 = as.numeric(prd_after['PVL_Immature_s1',]),
                        PVL = as.numeric(pvl_all),
                        PVLs3 = as.numeric(prd_after['PVL_Differentiated_s3',]))

# "T_cells_c3_CD4__Tfh_CXCL13"   T_cells_c8_CD8__LAG3
# "PVL_Differentiated_s3"          "PVL_Immature_s1"  "PVL_Immature_s2"

data_matrix[,'surstat']=as.character(Clinical[,'death'])
data_matrix[,'surtime']=as.numeric(Clinical[,'os'])
colnames(data_matrix)[5]='death'
# data_matrix[,'CX_PVL']<-data_matrix[,'PVLs1']/data_matrix[,'PVL']
data_matrix[,'CX_PVL']<-data_matrix[,'CXCL13']*data_matrix[,'PVLs1']/data_matrix[,'PVL']/data_matrix[,'PVLs3']

# /data_matrix[,'PVL']
data_matrix=data_matrix[!(data_matrix[,'death']=='0' & data_matrix[,'surtime']<730),]


data_matrix[,'datasets']='GSE58812'

saveRDS(data_matrix,file='~/ReCIDE/应用_TNBC/论文中结果/箱线图/GSE58812.rds')
