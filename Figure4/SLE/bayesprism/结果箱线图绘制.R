library(ggplot2)
library(ggpubr)

prd_df <- readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bayesprism/prd_df.rds")
prd_df=prd_df[,sort(colnames(prd_df))]

metadata<- readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bulkdata/GSE50772_meta.rds")
metadata=metadata[sort(row.names(metadata)),]
# [1] "B_exhausted"           "B_immature"           
# [3] "B_naive"               "B_non_switched_memory"
# [5] "B_switched_memory"     "CD14_mono"            
# [7] "CD16_mono"             "CD4_CM"               
# [9] "CD4_EM"                "CD4_IL22"             
# [11] "CD4_Naive"             "CD4_Th"               
# [13] "CD8_EM"                "CD8_Naive"            
# [15] "CD8_TE"                "DC"                   
# [17] "HSC"                   "ILC"                  
# [19] "Lymph_prolif"          "MAIT"                 
# [21] "NK"                    "NKT"                  
# [23] "Plasma"                "Platelets"            
# [25] "RBC"                   "Treg"                 
# [27] "gdT"                   "pDC"                  
# [29] "category"             

all(row.names(metadata)==colnames(prd_df))

plot_data=as.data.frame(t(prd_df))

plot_data[,'category']=metadata[,'disease.status.ch1']

ggplot(plot_data, aes(x= category, y=MAIT,fill=category)) + 
  geom_boxplot(aes(color = category),#这里的fill如果不设就是空心的
               size = .7, alpha = .5, outlier.size = 1,
               position = position_dodge(width = 0.5), width = 0.5) +
  # position_dodge是箱子间距离  width是箱宽
  # geom_boxplot(data=subset(plot_data, methods %in% 'empty'), width = .5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, face = "bold", angle = -45),
    axis.text.y = element_text(size = 9, face = "bold"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    legend.position = 'none'
    # axis.title = element_text(size = 8)
  )+
  stat_compare_means()+
  geom_hline(yintercept = 0.001,col = 'Black',linewidth=0.6,linetype=2)+
  geom_hline(yintercept = 0.013,col = 'Black',linewidth=0.6,linetype=2)



