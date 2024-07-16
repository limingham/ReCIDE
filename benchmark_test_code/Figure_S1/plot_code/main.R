
library(dplyr)
library(ggplot2)
library(ggpubr)

deone_PAN <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/deone_PAN.rds")
deone_kidney <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/deone_kidney.rds")
cross_PFC <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/cross_PFC.rds")
cross_PBMC <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/cross_PBMC.rds")
mult_PBMC <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/mult_PBMC.rds")
mult_ER <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/mult_ER.rds")
high_res_TNBC <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/high_res_TNBC.rds")
high_res_CRC <- readRDS("~/ReCIDE/benchmark测试_最终/绘制Figure1_优化前后比较_PCC/plotdata/high_res_CRC.rds")

deone_PAN[,'datasets']='leave_one_PAN'

deone_kidney[,'datasets']='leave_one_kidney'

cross_PFC[,'datasets']='cross_PFC'

cross_PBMC[,'datasets']='cross_PBMC'

mult_PBMC[,'datasets']='mult_PBMC'

mult_ER[,'datasets']='mult_ER'

high_res_TNBC[,'datasets']='high_res_TNBC'

high_res_CRC[,'datasets']='high_res_CRC'


plot_data=rbind(deone_PAN,
                deone_kidney,
                cross_PFC,
                cross_PBMC,
                mult_PBMC,
                mult_ER,
                high_res_TNBC,
                high_res_CRC)

# saveRDS(plot_data,file = '~/SWORD/作图/散点图_1/plot_data.rds')
# 
# 
# plot_data=readRDS('~/SWORD/作图/散点图_1/plot_RMSE.rds')

plot_data[,'identity']=paste0(plot_data[,'methods'],plot_data[,'datasets'])

plot_data[,'category']=rep(c('com','sep'),24)

sub1=subset(plot_data,category=='com')
sub2=subset(plot_data,category=='sep')

plot_all=sub1
plot_all[,'RMSE2']=sub2[,'RMSE']
# 
# minus=(plot_all[,'RMSE2']-plot_all[,'RMSE'])/plot_all[,'RMSE']
# cor(minus,plot_all[,'RMSE'])
# 
# plot_all[,'minus']=-minus
# plot_all[,'methods']=c(rep('CIBERSORT',4),rep('DWLS',4),rep('FARDEEP',4))

shape_level <- nlevels(as.factor(plot_all[["datasets"]]))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}

plot_all[,'methods']=factor(plot_all[,'methods'],levels=c('DWLS_com',
                                                          'CIBERSORT_com ',
                                                          'FARDEEP_com'))

plot_all[,'RMSE_enhance']=1-plot_all[,'RMSE2']/plot_all[,'RMSE']

group_by(plot_all,methods) %>% summarize_each(mean)


###5 7.5
ggplot(plot_all, aes(RMSE, RMSE2)) +
  # ggtitle("SWORD-DWLS/BayesPrism-Major")+
  geom_point(aes(fill = methods, color = methods, shape=datasets),size = 4, alpha = 5,stroke = 0.9)+
  # geom_point(aes(shape=Year),size=3)
  # scale_shape_manual(values=shapes)+
  scale_y_continuous(limits = c(0.4,1))+
  scale_x_continuous(limits = c(0.4,1))+
  scale_shape_manual(values=shapes)+
  # scale_y_continuous(limits = c(0.05,0.1))+
  # scale_x_continuous(limits = c(0.15,0.18))+
  geom_abline(intercept = 0, slope = 1, col = 'black',linewidth=0.7)+
  # stat_cor(data = scatter_data_Rare_Bayes[,c(1,2)],label.x = 0.05,label.y = 0.84, size=3)+
  theme_classic(base_size = 18)+
  theme(axis.title.x = element_text(family ='Times',size=11,color = 'black', face = "bold"),
        axis.title.y = element_text(family ='Times',size=11,color = 'black', face = "bold"),
        # axis.line = element_line(size = 1),
        axis.text.x = element_text(family ='Times',size=11,color = 'black', face = "bold"),
        axis.text.y = element_text(family ='Times',size=11,color = 'black', face = "bold"),
        axis.ticks.x = element_blank(),
        # axis.ticks = element_line(size = 1,color='black'),
        # plot.margin = margin(2,1,1,2,'cm'),
        plot.title = element_text(family ='Times',size=11,color = 'black', face = "bold",hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")
        # panel.border = element_rect(color = "black", size = 2, fill = NA)
  )+
  theme(
    legend.title= element_text(family ='Times',size=9,color = 'black', face = "bold"),
    legend.text = element_text(family ='Times',size=9,color = 'black', face = "bold")
  )






bar_data_init=group_by(plot_all,methods) %>% summarize_each(mean)

bar_data_init=as.data.frame(bar_data_init)

bar_data=data.frame(Value=c(bar_data_init[,'RMSE'],bar_data_init[,'RMSE2']),
                    methods=c('DWLS','CIBERSORT','FARDEEP','DWLS','CIBERSORT','FARDEEP'),
                    category=c('before','before','before','after','after','after'))


bar_data[,'methods']=factor(bar_data[,'methods'],levels=c('DWLS','FARDEEP','CIBERSORT'))


ggplot(bar_data,aes(x=methods,y=Value,fill=category))+
  geom_bar(stat = 'identity', 
           #柱状图位置并排:
           position = 'dodge', #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
           width = 0.8,      #设置柱子宽度,使变量之间分开
           color='white')+        
  scale_fill_manual(values = c('#ff6852','#56b4e9'))+
  # scale_color_manual(values = c('#ff6852','#56b4e9'))+
  theme_classic(base_size = 18)+
  theme(axis.title.x = element_text(family ='Arial',size=11,color = 'black', face = "bold"),
        axis.title.y = element_text(family ='Arial',size=11,color = 'black', face = "bold"),
        # axis.line = element_line(size = 1),
        axis.text.x = element_text(family ='Arial',size=11,color = 'black', face = "bold"),
        axis.text.y = element_text(family ='Arial',size=11,color = 'black', face = "bold"),
        axis.ticks.x = element_blank(),
        # axis.ticks = element_line(size = 1,color='black'),
        # plot.margin = margin(2,1,1,2,'cm'),
        plot.title = element_text(family ='Arial',size=11,color = 'black', face = "bold",hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")
        # panel.border = element_rect(color = "black", size = 2, fill = NA)
  )+
  theme(
    legend.title= element_text(family ='Arial',size=9,color = 'black', face = "bold"),
    legend.text = element_text(family ='Arial',size=9,color = 'black', face = "bold")
  )

##500 750
