library(ggpubr)
library(ggplot2)


prd_bisque <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bisque/prd_df.rds")
prd_bisque=as.data.frame(t(prd_bisque))
prd_bisque[,'DCs']=prd_bisque[,'DC']+prd_bisque[,'pDC']

# prd_music <- as.data.frame(t(readRDS("~/ReCIDE/应用_前二/COVID19/music/MuSiC_prd_df.rds")))
prd_music <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/music/prd_df.rds")
prd_music=as.data.frame(t(prd_music))
prd_music[,'DCs']=prd_music[,'DC']+prd_music[,'pDC']

prd_BayesPrism <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bayesprism/prd_df.rds")
prd_BayesPrism=as.data.frame(t(prd_BayesPrism))
prd_BayesPrism[,'DCs']=prd_BayesPrism[,'DC']+prd_BayesPrism[,'pDC']

prd_DWLS_sep <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/ReCIDE_DWLS/结果/prd_df.rds")
prd_DWLS_sep=as.data.frame(t(prd_DWLS_sep))
prd_DWLS_sep[,'DCs']=prd_DWLS_sep[,'DC']+prd_DWLS_sep[,'pDC']

prd_DWLS_com <- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/DWLS/prd_df.rds")
prd_DWLS_com=as.data.frame(t(prd_DWLS_com))
prd_DWLS_com[,'DCs']=prd_DWLS_com[,'DC']+prd_DWLS_com[,'pDC']

prd_SCDC<- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/scdc/prd_df.rds")
prd_SCDC=as.data.frame(t(prd_SCDC))
prd_SCDC[,'DCs']=prd_SCDC[,'DC']+prd_SCDC[,'pDC']
# prd_DWLS_com=subset(prd_DWLS_com,DC>0)
# prd_BayesPrism=subset(prd_BayesPrism,DC>0)

metadata<- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/metadata_query.rds")
metadata=metadata[sort(row.names(metadata)),]
patient_names=row.names(metadata)

prd_bisque <-prd_bisque[patient_names,sort(colnames(prd_bisque))]
prd_music <-prd_music[patient_names,sort(colnames(prd_music))]
prd_BayesPrism <-prd_BayesPrism[patient_names,sort(colnames(prd_BayesPrism))]
prd_DWLS_sep <-prd_DWLS_sep[patient_names,sort(colnames(prd_DWLS_sep))]
prd_DWLS_com <-prd_DWLS_com[patient_names,sort(colnames(prd_DWLS_com))]
prd_SCDC <-prd_SCDC[patient_names,sort(colnames(prd_SCDC))]

colnames(prd_DWLS_sep)=colnames(prd_bisque)
colnames(prd_DWLS_com)=colnames(prd_bisque)


prd_bisque[,'category']=metadata[,'characteristics_ch1.3']
prd_music[,'category']=metadata[,'characteristics_ch1.3']
prd_BayesPrism[,'category']=metadata[,'characteristics_ch1.3']
prd_DWLS_sep[,'category']=metadata[,'characteristics_ch1.3']
prd_DWLS_com[,'category']=metadata[,'characteristics_ch1.3']
prd_SCDC[,'category']=metadata[,'characteristics_ch1.3']


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

plot_data=subset(plot_data,category %in% c("sample timing: healthy",
                                           "sample timing: acute"))
# for (i in 1:nrow(plot_data)) {
#   if(plot_data[i,'category']=="status: CTL"){plot_data[i,'category']<-'Control'}
#   if(plot_data[i,'category']=="status: MCI"){plot_data[i,'category']<-'MCI'}
#   if(plot_data[i,'category']=="status: AD"){plot_data[i,'category']<-'AD'}
#   
# }
# # 
# 
plot_data[,'category']=factor(plot_data[,'category'],levels=c("sample timing: healthy","sample timing: acute"))
plot_data[,'methods']=factor(plot_data[,'methods'],levels=c("ReCIDE-DWLS","DWLS","BayesPrism","Bisque","MuSiC",'SCDC'))

library(RColorBrewer)
vec_color <- c('#a52a2a','#4861a3')


# stat.test <- plot_data %>%
#   group_by(methods) %>%
#   pairwise_wilcox_test(
#     DC ~ category, paired = FALSE,
#     p.adjust.method = "fdr"
#   )
# wilcox.test(stat.test[,'B'] ~ stat.test[,'category'])
# pairwise_wilcox_test
# stat.test <- stat.test %>% add_xy_position(x = "methods")
# stat.test=as.data.frame(stat.test)
# colnames(stat.test)[12]='category'

ggplot(plot_data, aes(x=methods, y=DCs, fill=category,color = category)) + 
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
  stat_compare_means(label.y=0.07,method = "wilcox.test")+
  geom_hline(yintercept = 0.005,col = 'Black',linewidth=0.6,linetype=2)+
  geom_hline(yintercept = 0.02,col = 'Black',linewidth=0.6,linetype=2)

