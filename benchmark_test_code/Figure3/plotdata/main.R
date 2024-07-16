library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####CRC
#######

prd_BayesPrism<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/bayesprism/prd_df.rds")
prd_Bisque<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/bisque/prd_df.rds")

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/DWLS_com/prd_df.rds")
prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/DWLS_sep/prd_df.rds")

prd_MuSiC<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/music/prd_df.rds")
prd_SCDC <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/scdc/prd_df.rds")


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

key_df <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_CRC/EXP_and_KEY/KEY.rds")
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]


all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('CRC_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_high_res_CRC.rds')



#######
####TNBC
#######
# prd_bisque_df <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/bisque/prd_bisque_df.rds")
# prd_bisque_df['NEU',]=0
# saveRDS(prd_bisque_df,file="~/ReCIDE/benchmark测试_最终/deone_kidney/bisque/prd_bisque_df.rds")
key_df <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/EXP_and_KEY/KEY.rds")
row.names(key_df) <- str_replace(row.names(key_df), "\\+", "_")
row.names(key_df) <- str_replace(row.names(key_df), "-", "_")
row.names(key_df) <- str_replace(row.names(key_df), " ", "_")
row.names(key_df) <- str_replace(row.names(key_df), " ", "_")
row.names(key_df) <- str_replace(row.names(key_df), " ", "_")

key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/DWLS_com/prd_df.rds")
prd_DWLS_com[is.na(prd_DWLS_com)]=0
# 
# prd_DWLS_com_more=row.names(prd_DWLS_com)[!(row.names(prd_DWLS_com) %in% row.names(key_df))]
# key_df_more=row.names(key_df)[!(row.names(key_df) %in% row.names(prd_DWLS_com))]
# 
# prd_DWLS_com[key_df_more,]=0
# key_df[prd_DWLS_com_more,]=0


prd_BayesPrism <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/prd_df.rds")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), "\\+", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), "-", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), " ", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), " ", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), " ", "_")
prd_BayesPrism[is.na(prd_BayesPrism)]=0
# prd_BayesPrism[prd_DWLS_com_more,]=0


prd_Bisque<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bisque/prd_df.rds")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), "\\+", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), "-", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), " ", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), " ", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), " ", "_")
prd_Bisque[is.na(prd_Bisque)]=0
# prd_Bisque[prd_DWLS_com_more,]=0



prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/DWLS_sep/prd_df.rds")
# prd_DWLS_sep[prd_DWLS_com_more,]=0

prd_MuSiC  <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/music/prd_df.rds")
prd_MuSiC['Myeloid_c5_Macrophage_3_SIGLEC1',]=0
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), "\\+", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), "-", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), " ", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), " ", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), " ", "_")
prd_MuSiC[is.na(prd_MuSiC)]=0
# prd_MuSiC[prd_DWLS_com_more,]=0

prd_SCDC <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/scdc/prd_df.rds")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), "\\+", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), "-", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), " ", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), " ", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), " ", "_")
prd_SCDC[is.na(prd_SCDC)]=0
# prd_SCDC[prd_DWLS_com_more,]=0


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]



all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('TNBC_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_high_res_TNBC.rds')

library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####PAN
#######

prd_BayesPrism<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/bayesprism/prd_df.rds")
prd_Bisque <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/bisque/prd_df.rds")

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/DWLS_com/prd_df.rds")
prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/DWLS_sep/prd_df.rds")

prd_MuSiC <- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/music/prd_df.rds")
prd_SCDC<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/scdc/prd_df.rds")

prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

key_df<- readRDS("~/ReCIDE/benchmark测试_最终/deone_PAN/EXP_and_KEY/KEY.rds")
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]


all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('PAN_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_high_res_PAN.rds')




library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####kidney
#######

prd_BayesPrism <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/bayesprism/prd_df.rds")
prd_BayesPrism['NEU',]=0
prd_Bisque <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/bisque/prd_df.rds")

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/DWLS_com/prd_df.rds")
prd_DWLS_com['NEU',]=0

prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/DWLS_sep/prd_df.rds")
prd_DWLS_sep['NEU',]=0

prd_MuSiC <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/music/prd_df.rds")
prd_SCDC <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/scdc/prd_df.rds")


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

key_df <- readRDS("~/ReCIDE/benchmark测试_最终/deone_kidney/EXP_and_KEY/KEY_kidney.rds")
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]


all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('Kidney_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_high_res_Kidney.rds')




library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####mult_pbmc
#######

prd_BayesPrism<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/bayesprism/prd_df.rds")
prd_Bisque <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/bisque/prd_df.rds")

prd_DWLS_com<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/DWLS_com/prd_df.rds")
prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/DWLS_sep/prd_df.rds")

prd_MuSiC <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/music/prd_df.rds")
prd_SCDC<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/scdc/prd_df.rds")


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

key_df<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_PBMC/EXP_and_KEY/key_PBMC.rds")
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]


all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('mult_pbmc_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_mult_pbmc.rds')




library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####mult_ER
#######
prd_BayesPrism<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/bayesprism/prd_df.rds")
prd_Bisque <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/bisque/prd_df.rds")

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/DWLS_com/prd_df.rds")
prd_DWLS_sep<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/DWLS_sep/prd_df.rds")

prd_MuSiC <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/music/prd_df.rds")
prd_SCDC <- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/scdc/prd_df.rds")


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

key_df<- readRDS("~/ReCIDE/benchmark测试_最终/mult_reference_ER新/EXP_and_KEY/KEY.rds")
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]



all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('mult_ER_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_mult_ER.rds')




library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####cross_pfc
#######

prd_BayesPrism <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/bayesprism/prd_df.rds")
prd_Bisque<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/bisque/prd_df.rds")

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/DWLS_com/prd_df.rds")
prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/DWLS_sep/prd_df.rds")

prd_MuSiC <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/music/prd_df.rds")
prd_SCDC  <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/scdc/prd_df.rds")


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),sort(colnames(prd_Bisque))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),sort(colnames(prd_MuSiC))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

key_df<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PFC/EXP_and_key/key_PFC.rds")
key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]


all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('cross_pfc_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_cross_pfc.rds')




library(ggplot2)
library(ggpubr)
library(ggpubr)
library(dplyr)
library(reshape2)
#######
####cross_pbmc
#######
patient_names=c("CV0902","CV0904","CV0911","CV0915","CV0917","CV0926","CV0929",     
                "CV0934","CV0939","CV0940","CV0944","MH8919176","MH8919177","MH8919178",  
                "MH8919179","MH8919226","MH8919227","MH8919282","MH8919283","MH8919332","MH8919333",  
                "newcastle65","newcastle74")

prd_BayesPrism <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/BayesPrism/prd_df.rds")
prd_Bisque<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/bisque/prd_df.rds")

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/DWLS_com/prd_df.rds")
prd_DWLS_sep<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/DWLS_sep/prd_df.rds")

prd_MuSiC<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/music/prd_df.rds")
prd_SCDC<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/scdc/prd_df.rds")


prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),patient_names]
prd_Bisque=prd_Bisque[sort(row.names(prd_Bisque)),patient_names]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),patient_names]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),patient_names]
prd_MuSiC=prd_MuSiC[sort(row.names(prd_MuSiC)),patient_names]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),patient_names]

key_df<- readRDS("~/ReCIDE/benchmark测试_最终/cross_dataset_PBMC/EXP_and_KEY/key_combine.rds")
key_df=key_df[sort(row.names(key_df)),patient_names]


all(row.names(prd_BayesPrism)==row.names(key_df))
all(row.names(prd_Bisque)==row.names(key_df))
all(row.names(prd_DWLS_com)==row.names(key_df))
all(row.names(prd_DWLS_sep)==row.names(key_df))
all(row.names(prd_MuSiC)==row.names(key_df))
all(row.names(prd_SCDC)==row.names(key_df))

all(colnames(prd_BayesPrism)==colnames(key_df))
all(colnames(prd_Bisque)==colnames(key_df))
all(colnames(prd_DWLS_com)==colnames(key_df))
all(colnames(prd_DWLS_sep)==colnames(key_df))
all(colnames(prd_MuSiC)==colnames(key_df))
all(colnames(prd_SCDC)==colnames(key_df))


prd_BayesPrism_melt=melt(prd_BayesPrism)
prd_BayesPrism_melt[,'celltype']=rep(row.names(prd_BayesPrism),ncol(prd_BayesPrism))
prd_BayesPrism_melt[,'true_prop']=melt(key_df)[,2]
prd_BayesPrism_melt[,'methods']='BayesPrism'

prd_Bisque_melt=melt(prd_Bisque)
prd_Bisque_melt[,'celltype']=rep(row.names(prd_Bisque),ncol(prd_Bisque))
prd_Bisque_melt[,'true_prop']=melt(key_df)[,2]
prd_Bisque_melt[,'methods']='Bisque'

prd_DWLS_com_melt=melt(prd_DWLS_com)
prd_DWLS_com_melt[,'celltype']=rep(row.names(prd_DWLS_com),ncol(prd_DWLS_com))
prd_DWLS_com_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_com_melt[,'methods']='DWLS'


prd_DWLS_sep_melt=melt(prd_DWLS_sep)
prd_DWLS_sep_melt[,'celltype']=rep(row.names(prd_DWLS_sep),ncol(prd_DWLS_sep))
prd_DWLS_sep_melt[,'true_prop']=melt(key_df)[,2]
prd_DWLS_sep_melt[,'methods']='ReCIDE-DWLS'

prd_MuSiC_melt=melt(prd_MuSiC)
prd_MuSiC_melt[,'celltype']=rep(row.names(prd_MuSiC),ncol(prd_MuSiC))
prd_MuSiC_melt[,'true_prop']=melt(key_df)[,2]
prd_MuSiC_melt[,'methods']='MuSiC'

prd_SCDC_melt=melt(prd_SCDC)
prd_SCDC_melt[,'celltype']=rep(row.names(prd_SCDC),ncol(prd_SCDC))
prd_SCDC_melt[,'true_prop']=melt(key_df)[,2]
prd_SCDC_melt[,'methods']='SCDC'


plot_data=rbind(prd_BayesPrism_melt,
                prd_Bisque_melt,
                prd_DWLS_com_melt,
                prd_DWLS_sep_melt,
                prd_MuSiC_melt,
                prd_SCDC_melt)



colnames(plot_data)=c('patient','prd_prop','celltype','true_prop','methods')
plot_data[,'category']=paste0('cross_pbmc_',plot_data[,'celltype'])
saveRDS(plot_data,file = '~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/plot_data_cross_pbmc.rds')
