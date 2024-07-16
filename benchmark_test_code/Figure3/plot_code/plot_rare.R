library(dplyr)
library(ggpubr)
library(RColorBrewer)


# benchmark测试_最终
file_list=list.files('~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata')
file_list=file_list[file_list %in% c(
  "plot_data_cross_pbmc.rds","plot_data_cross_pfc.rds",
  "plot_data_deone_Kidney.rds","plot_data_deone_PAN.rds","plot_data_high_res_CRC.rds",
  "plot_data_high_res_TNBC.rds", "plot_data_mult_ER.rds","plot_data_mult_pbmc.rds")]

# file_list=file_list[file_list %in% c(
#   "plot_data_high_res_CRC.rds","plot_data_high_res_TNBC.rds", "plot_data_mult_ER.rds")]
# 

combine_data=readRDS(paste0('~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/',file_list[1]))
for (i in 2:length(file_list)) {
  data_in=readRDS(paste0('~/ReCIDE/benchmark测试_最终/绘制Figure1_回归散点图/plotdata/',file_list[i]))
  combine_data=rbind(combine_data,data_in)
}

# saveRDS(combine_data,file = '~/SWORD/作图/稀有评估散点图/plot_code/处理好数据/combine_data.rds')

combine_data[,'ident']=combine_data[,'category']
data_mean = group_by(combine_data,ident) %>% summarize_each(mean)
data_mean=as.data.frame(data_mean)
row.names(data_mean)=data_mean[,1]

# data_mean_filter=data_mean[data_mean[,'true_prop']<0.1 & data_mean[,'true_prop']>0.02,]
# plotdata=subset(plotdata,prd_prop<0.25 & true_prop<0.25)
data_mean_filter=data_mean[data_mean[,'true_prop']<0.02,]

plotdata=subset(combine_data,ident %in% row.names(data_mean_filter))
table(row.names(data_mean_filter))

#####线性回归并绘制柱状图
#####线性回归并绘制柱状图
slope_c=c(
  lm(plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS','true_prop'] ~ plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS','prd_prop'] + 0)[["coefficients"]][1],
  lm(plotdata[plotdata[,'methods'] == 'SCDC','true_prop'] ~ plotdata[plotdata[,'methods'] == 'SCDC','prd_prop'] + 0)[["coefficients"]][1],
  lm(plotdata[plotdata[,'methods'] == 'DWLS','true_prop'] ~ plotdata[plotdata[,'methods'] == 'DWLS','prd_prop'] + 0)[["coefficients"]][1],
  lm(plotdata[plotdata[,'methods'] == 'BayesPrism','true_prop'] ~ plotdata[plotdata[,'methods'] == 'BayesPrism','prd_prop'] + 0)[["coefficients"]][1],
  lm(plotdata[plotdata[,'methods'] == 'MuSiC','true_prop'] ~ plotdata[plotdata[,'methods'] == 'MuSiC','prd_prop'] + 0)[["coefficients"]][1],
  lm(plotdata[plotdata[,'methods'] == 'Bisque','true_prop'] ~ plotdata[plotdata[,'methods'] == 'Bisque','prd_prop'] + 0)[["coefficients"]][1]
)
# slope_c=c(ModelMetrics::rmse(plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS','true_prop'] , plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS','prd_prop']),
#           ModelMetrics::rmse(plotdata[plotdata[,'methods'] == 'DWLS','true_prop'] , plotdata[plotdata[,'methods'] == 'DWLS','prd_prop']),
#           ModelMetrics::rmse(plotdata[plotdata[,'methods'] == 'BayesPrism','true_prop'] , plotdata[plotdata[,'methods'] == 'BayesPrism','prd_prop']),
#           ModelMetrics::rmse(plotdata[plotdata[,'methods'] == 'Bisque','true_prop'] , plotdata[plotdata[,'methods'] == 'Bisque','prd_prop'] ),
#           ModelMetrics::rmse(plotdata[plotdata[,'methods'] == 'MuSiC','true_prop'] , plotdata[plotdata[,'methods'] == 'MuSiC','prd_prop'] ),
#           ModelMetrics::rmse(plotdata[plotdata[,'methods'] == 'SCDC','true_prop'] , plotdata[plotdata[,'methods'] == 'SCDC','prd_prop'] ))

# slope_c=c(cor(plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS','true_prop'] , plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS','prd_prop']),
#           cor(plotdata[plotdata[,'methods'] == 'DWLS','true_prop'] , plotdata[plotdata[,'methods'] == 'DWLS','prd_prop']),
#           cor(plotdata[plotdata[,'methods'] == 'BayesPrism','true_prop'] , plotdata[plotdata[,'methods'] == 'BayesPrism','prd_prop']),
#           cor(plotdata[plotdata[,'methods'] == 'Bisque','true_prop'] , plotdata[plotdata[,'methods'] == 'Bisque','prd_prop'] ),
#           cor(plotdata[plotdata[,'methods'] == 'MuSiC','true_prop'] , plotdata[plotdata[,'methods'] == 'MuSiC','prd_prop'] ),
#           cor(plotdata[plotdata[,'methods'] == 'SCDC','true_prop'] , plotdata[plotdata[,'methods'] == 'SCDC','prd_prop'] ))



slope_df<-as.data.frame(slope_c)
slope_df[,2]<-c('ReCIDE-DWLS','SCDC','DWLS','BayesPrism','MuSiC','Bisque')

vec.color <-rev(c('#a52a2a','#56b4e9','#F58840','#008000','#5b4e59',"#761c77"))

# slope_df[,1]=abs(slope_df[,1]-1)
colnames(slope_df)<-c('value','methods')

# slope_df=slope_df %>% arrange(value)
slope_df[,'methods']=factor(slope_df[,'methods'],levels=slope_df[,'methods'])
###3，4.5
# plot.bar <- 
p1=ggplot(slope_df, aes(x = reorder(methods,value), y = value, color = methods, fill = methods)) +
  geom_bar(stat = 'identity') +
  geom_abline(intercept = 1, slope = 0, col = 'Black',linewidth=1,linetype=2)+
  scale_color_manual(values = rev(vec.color), labels = slope_df[,'methods']) +
  scale_fill_manual(values = rev(vec.color), labels = slope_df[,'methods']) +
  labs(title = "", y = 'Median of ranks', x = '') + 
  coord_flip() +
  theme_bw() +
  theme(panel.background = element_rect(color = 'black', size = 1,
                                        fill = 'transparent'),
        panel.border = element_rect(size = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = 'lightgray',
                                          size = 0.15, linetype = 'dashed'),
        legend.position = 'none',
        axis.text.x = element_text(size = 11, color = 'black', family = 'Helvetica'),
        axis.text.y = element_text(size = 11, color = 'black', family = 'Helvetica'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, color = 'black', family = 'Helvetica'))




plotdata[,'methods']=factor(plotdata[,'methods'],levels=rev(slope_df[,'methods']))

#####线性回归并绘制散点图

p2=ggplot(plotdata, aes(prd_prop, true_prop)) +
  # ggtitle("DWLS-Rare")+
  # ggtitle("DWLS-Minor")+
  ggtitle("DWLS-Major")+
  geom_point(aes(fill = methods, color = methods, shape=methods),size = 0.7, alpha = 0.7)+
  geom_abline(intercept = 0, slope = 1, col = 'Black',linewidth=0.6,linetype=2)+
  scale_color_manual(values = vec.color) + 
  scale_fill_manual(values = vec.color) + 
  scale_shape_manual(values = c(1,2,3,4,5,6))+
  coord_cartesian(xlim = c(0.002, 0.06), ylim = c(0.0015, 0.06))+
  geom_smooth(data = plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[6],linewidth=0.6)+
  geom_smooth(data = plotdata[plotdata[,'methods'] == 'SCDC',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[5],linewidth=0.6)+
  geom_smooth(data = plotdata[plotdata[,'methods'] == 'DWLS',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[4],linewidth=0.6)+
  geom_smooth(data = plotdata[plotdata[,'methods'] == 'BayesPrism',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[3],linewidth=0.6)+
  geom_smooth(data = plotdata[plotdata[,'methods'] == 'MuSiC',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[2],linewidth=0.6)+
  geom_smooth(data = plotdata[plotdata[,'methods'] == 'Bisque',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[1],linewidth=0.6)+
  # geom_smooth(data = plotdata[plotdata[,'methods'] == 'ReCIDE-DWLS',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[1],linewidth=0.45,fill = vec.color[1])+
  # geom_smooth(data = plotdata[plotdata[,'methods'] == 'DWLS',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[2],linewidth=0.45,fill = vec.color[2])+
  # geom_smooth(data = plotdata[plotdata[,'methods'] == 'BayesPrism',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[3],linewidth=0.45,fill = vec.color[3])+
  # geom_smooth(data = plotdata[plotdata[,'methods'] == 'Bisque',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[4],linewidth=0.45,fill = vec.color[4])+
  # geom_smooth(data = plotdata[plotdata[,'methods'] == 'MuSiC',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[5],linewidth=0.45,fill = vec.color[5])+
  # geom_smooth(data = plotdata[plotdata[,'methods'] == 'SCDC',c('prd_prop','true_prop')],method = 'lm', formula = y ~ x + 0, se = FALSE, fullrange = TRUE,col = vec.color[6],linewidth=0.45,fill = vec.color[6])+
  # # coord_cartesian(xlim = c(0.0015, 0.05), ylim = c(0.0015, 0.05))
  theme_classic(base_size = 18)+
  theme(axis.title.x = element_text(family ='Helvetica',size=13,color = 'black'),
        axis.title.y = element_text(family ='Helvetica',size=13,color = 'black'),
        # axis.line = element_line(size = 1),
        axis.text.x = element_text(family ='Helvetica',size=13,color = 'black'),
        axis.text.y = element_text(family ='Helvetica',size=13,color = 'black'),
        axis.ticks.x = element_blank(),
        # axis.ticks = element_line(size = 1,color='black'),
        # plot.margin = margin(2,1,1,2,'cm'),
        # plot.title = element_blank(),
        plot.title = element_text(family ='Helvetica',size=11,color = 'black',hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")
        # panel.border = element_rect(color = "black", size = 2, fill = NA)
        # xlim=c(-0.004,0.1),ylim=c(-0.003,0.1)
  )+
  theme(
    legend.title= element_text(family ='Helvetica',size=10,color = 'black'),
    legend.text = element_text(family ='Helvetica',size=10,color = 'black'),

  )+
  guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1),
         fill = guide_legend(ncol = 1))+
theme(legend.key.size = unit(1, 'cm'))+ 
guides(color = guide_legend(override.aes = list(size = 4)))

p2+p1




##7.5 5.5 3 5.3