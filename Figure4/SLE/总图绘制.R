library(ggplot2)
library(ggpubr)

# prd_bisque <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/SLE_true_new/bisque/bisque_prd_df.rds")))
# prd_BayesPrism <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/SLE_true_new/BayesPrism/Bayes_prd_df.rds")))
# prd_DWLS_sep <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/SLE_true_new/ReCIDE_DWLS/prd_sep_LowToHigh_PCA.rds")))
# prd_DWLS_com <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/SLE_true_new/DWLS_build/prd_com_PBMC.rds")))
# prd_SCDC <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/SLE_true_new/scdc_原/SCDC_prd_df.rds")))
# prd_music <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/SLE_true_new/music_原/MuSiC_prd_df.rds")))
# 

prd_bisque <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bisque/prd_df.rds")))
prd_BayesPrism <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bayesprism/prd_df.rds")))
prd_DWLS_sep <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二_新_inter/SLE/ReCIDE_DWLS/结果/prd_df.rds")))
prd_DWLS_com <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二_新_inter/SLE/DWLS/prd_df.rds")))
prd_SCDC <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二_新_inter/SLE/scdc/prd_df.rds")))
prd_music <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二_新_inter/SLE/music/prd_df.rds")))


prd_bisque=prd_bisque[sort(row.names(prd_bisque)),sort(colnames(prd_bisque))]
prd_music=prd_music[sort(row.names(prd_music)),sort(colnames(prd_music))]
prd_BayesPrism=prd_BayesPrism[sort(row.names(prd_BayesPrism)),sort(colnames(prd_BayesPrism))]
prd_DWLS_sep=prd_DWLS_sep[sort(row.names(prd_DWLS_sep)),sort(colnames(prd_DWLS_sep))]
prd_DWLS_com=prd_DWLS_com[sort(row.names(prd_DWLS_com)),sort(colnames(prd_DWLS_com))]
prd_SCDC=prd_SCDC[sort(row.names(prd_SCDC)),sort(colnames(prd_SCDC))]

colnames(prd_DWLS_sep)=colnames(prd_BayesPrism)
colnames(prd_DWLS_com)=colnames(prd_BayesPrism)

# prd_BayesPrism=subset(prd_BayesPrism,Lymph.prolif>0)


metadata<- readRDS("~/ReCIDE/应用_前二_新_inter/SLE/bulkdata/GSE50772_meta.rds")
metadata=metadata[sort(row.names(metadata)),]

patient_names=row.names(metadata)

prd_bisque <-prd_bisque[patient_names,]
prd_music <-prd_music[patient_names,]
prd_BayesPrism <-prd_BayesPrism[patient_names,]
prd_DWLS_sep <-prd_DWLS_sep[patient_names,]
prd_DWLS_com <-prd_DWLS_com[patient_names,]
prd_SCDC <-prd_SCDC[patient_names,]


prd_bisque[,'category']=metadata[,'disease.status.ch1']
prd_music[,'category']=metadata[,'disease.status.ch1']
prd_BayesPrism[,'category']=metadata[,'disease.status.ch1']
prd_DWLS_sep[,'category']=metadata[,'disease.status.ch1']
prd_DWLS_com[,'category']=metadata[,'disease.status.ch1']
prd_SCDC[,'category']=metadata[,'disease.status.ch1']

# prd_DWLS_com=subset(prd_DWLS_com,Lymph.prolif>0)
# prd_music=subset(prd_music,Lymph.prolif>0)

prd_bisque[,'methods']='Bisque'
prd_music[,'methods']='MuSiC'
prd_BayesPrism[,'methods']='BayesPrism'
prd_DWLS_sep[,'methods']='ReCIDE-DWLS'
prd_DWLS_com[,'methods']='DWLS'
prd_SCDC[,'methods']='SCDC'

plot_data=rbind(prd_bisque,
                prd_music,
                prd_BayesPrism,
                prd_DWLS_sep,
                prd_DWLS_com,
                prd_SCDC)

# plot_data=subset(plot_data,Lymph.prolif>0)


# 
plot_data[,'category']=factor(plot_data[,'category'],levels=c("Control","SLE"))
plot_data[,'methods']=factor(plot_data[,'methods'],levels=c("ReCIDE-DWLS","DWLS","BayesPrism","Bisque","MuSiC","SCDC"))

library(RColorBrewer)
vec_color <- c('#a52a2a','#4861a3')

# library(rstatix)
# stat.test <- plot_data %>%
#   group_by(methods) %>%
#   pairwise_wilcox_test(
#     Lymph.prolif ~ category, paired = FALSE,
#     p.adjust.method = 'fdr'
#   )
# # wilcox.test(stat.test[,'Lymph.prolif'] ~ stat.test[,'category'])
# # pairwise_wilcox_test
# stat.test <- stat.test %>% add_xy_position(x = "methods")
# stat.test=as.data.frame(stat.test)
# colnames(stat.test)[12]='category'

ggplot(plot_data, aes(x=methods, y=MAIT, fill=category,color = category)) + 
  geom_boxplot(width = .6, size = .7, alpha = .5, outlier.size = 1, ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, face = "plain", angle = -45),
    axis.text.y = element_text(size = 10, face = "plain"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.4, "cm"),
    # axis.title = element_text(size = 8)
    legend.position = "bottom"
  )+
  labs (x = "", y = "") +
  coord_cartesian(ylim = c(0, 0.08))+
  scale_color_manual(values= vec_color) +
  # guides(fill = "none", color = "none")+
  scale_fill_manual(values= vec_color)+
  # geom_hline(yintercept = 0.02,col = 'Black',linewidth=0.6,linetype=2)+
  geom_hline(yintercept = 0.045,col = 'Black',linewidth=0.6,linetype=2)+
  geom_hline(yintercept = 0.00015,col = 'Black',linewidth=0.6,linetype=2)+
  # geom_hline(yintercept = 0.15,col = 'Black',linewidth=0.6,linetype=2)+
  stat_compare_means(label.y=0.06,method = "wilcox.test")