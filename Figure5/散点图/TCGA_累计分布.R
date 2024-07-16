library(survival)
library(survminer)
prd_after<- readRDS("~/SWORD/应用/TNBC/TCGA/SWORD_DWLS_results/prd_sep_df_ref150_30.rds")
load("~/SWORD/应用/TNBC/TCGA/target_bulk/Clinical_all.RData")
Clinical=phonedata[colnames(prd_after),]

pvl_all=prd_after['PVL_Differentiated_s3',]+prd_after['PVL_Immature_s1',]+prd_after['PVL_Immature_s2',]+prd_after['Cycling_PVL',]

data_matrix<-data.frame(CXCL13 = as.numeric(prd_after['T_cells_c3_CD4__Tfh_CXCL13',]),
                        PVLs1 = as.numeric(prd_after['PVL_Immature_s1',]),
                        PVL = as.numeric(pvl_all),
                        PVLs3 = as.numeric(prd_after['PVL_Differentiated_s3',]))


# "T_cells_c3_CD4__Tfh_CXCL13"   T_cells_c8_CD8__LAG3
# "PVL_Differentiated_s3"          "PVL_Immature_s1"  "PVL_Immature_s2"

data_matrix[,'surstat']=as.character(Clinical[,'death'])
data_matrix[,'surtime']=as.numeric(Clinical[,'os'])


# plot_scatter<-subset(data_matrix,(surstat == 0 & surtime<1825))
plot_scatter<-data_matrix
plot_scatter[,'ecdf']=plot_scatter[,'CXCL13']/as.numeric(pvl_all)
# /plot_scatter[,'PVL']
plot_scatter[,'category']='live'

plot_scatter[plot_scatter[,'surstat']=="1",'category']='death'

plot_scatter=plot_scatter[!(plot_scatter[,'category']=='live' & plot_scatter[,'surtime']<730),]





ggplot(plot_scatter,
       aes(
         x=ecdf,
         color=category
       ))+
  stat_ecdf( # ggplot2中的经验累积分布函数
    size=1   # 线条粗细
  )+
  theme_classic(base_size = 18)+
  scale_color_manual(values = c("#348c4a", "#9d4649"))+
  theme(axis.title.x = element_text(family ='Helvetica',size=13,color = 'black'),
        axis.title.y = element_text(family ='Helvetica',size=13,color = 'black'),
        # axis.line = element_line(size = 1),
        axis.text.x = element_text(family ='Helvetica',size=13,color = 'black'),
        axis.text.y = element_text(family ='Helvetica',size=13,color = 'black'),
        axis.ticks.x = element_blank(),
        # axis.ticks = element_line(size = 1,color='black'),
        # plot.margin = margin(2,1,1,2,'cm'),
        # plot.title = element_blank(),
        plot.title = element_text(family ='Helvetica',size=13,color = 'black',hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")
        # panel.border = element_rect(color = "black", size = 2, fill = NA)
        # xlim=c(-0.004,0.1),ylim=c(-0.003,0.1)
  )+
  theme(
    legend.title= element_text(family ='Helvetica',size=9,color = 'black'),
    legend.text = element_text(family ='Helvetica',size=9,color = 'black'),
    # legend.key.size = unit(35, "pt"),
    # legend.key.height = unit(35, "pt"),
    # legend.key.width = unit(55, "pt")
  )+
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1),
         fill = guide_legend(ncol = 1))+
  theme(legend.key.size = unit(1, 'cm'))+ 
  guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_vline(xintercept=0.043,linetype=3)
  
####5.8   3
