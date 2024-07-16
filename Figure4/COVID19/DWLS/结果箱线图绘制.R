library(ggplot2)
library(ggpubr)

prd_df<- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/DWLS/prd_com_PBMC.rds")
# prd_df['B',]=prd_df[2,]+prd_df[3,]+prd_df[4,]
prd_df=prd_df[,sort(colnames(prd_df))]
prd_df['DCs',]=prd_df['DC',]+prd_df['pDC',]

metadata<- readRDS("~/ReCIDE/应用_前二_新_inter/COVID19/bulkdata/metadata_query.rds")
metadata=metadata[sort(row.names(metadata)),]


all(row.names(metadata)==colnames(prd_df))

plot_data=as.data.frame(t(prd_df))

plot_data[,'category']=metadata[,'characteristics_ch1.3']
# 
# patient_names=row.names(metadata[metadata[,"included in case -control study:ch1"]=='yes',])
# plot_data=plot_data[patient_names,]

plot_data=subset(plot_data,category %in% c("sample timing: healthy",
                                           "sample timing: acute"))
# # # 
# # 
# plot_data[,'category']=factor(plot_data[,'category'],levels = c("sample timing: healthy",
#                                                                 "sample timing: acute"))


ggplot(plot_data, aes(x= category, y=DCs,fill=category)) + 
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
  stat_compare_means()
  
colnames(plot_data)
cor_ref=c(0.1228707508,0.0787832509,0.0010085376,0.0006977433,0.0038242098,0.3803133465,
  0.1825432436,0.0154265085,0.1773708866,0.0310588410,0.0061026815)

cor_query=plot_data[plot_data[,'category'] == "status: CTL",1:11]



cor_resuts=c()
for (i in 1:nrow(cor_query)) {
  cor_resuts[i]=cor(cor_ref,as.numeric(cor_query[i,]))
  
}




