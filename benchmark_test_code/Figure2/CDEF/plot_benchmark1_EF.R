library(dplyr)
library(ggplot2)
library(ggpubr)

source("~/ReCIDE/benchmark测试/绘制Figure1_箱线图/mult_ER新/main.R")
plot_data1=subset(plot_data,methods %in% c('BayesPrism',
                                           "Bisque",
                                           "DWLS_sep",
                                           "DWLS_com",
                                           "MuSiC",
                                           "SCDC"))
plot_data1[,'datasets']='deone_kidney'


source("~/ReCIDE/benchmark测试/绘制Figure1_箱线图/mult_PBMC/main.R")
plot_data2=subset(plot_data,methods %in% c('BayesPrism',
                                           "Bisque",
                                           "DWLS_sep",
                                           "DWLS_com",
                                           "MuSiC",
                                           "SCDC"))
plot_data2[,'datasets']='deone_PFC'


source("~/ReCIDE/benchmark测试/绘制Figure1_箱线图/high_res_CRC/main.R")
plot_data3=subset(plot_data,methods %in% c('BayesPrism',
                                           "Bisque",
                                           "DWLS_sep",
                                           "DWLS_com",
                                           "MuSiC",
                                           "SCDC"))
plot_data3[,'datasets']='cross_PFC'


source("~/ReCIDE/benchmark测试/绘制Figure1_箱线图/high_res_TNBC/main.R")
plot_data4=subset(plot_data,methods %in% c('BayesPrism',
                                           "Bisque",
                                           "DWLS_sep",
                                           "DWLS_com",
                                           "MuSiC",
                                           "SCDC"))
plot_data4[,'datasets']='cross_PBMC'

plot_data=rbind(plot_data1,
                plot_data2,
                plot_data3,
                plot_data4)

# saveRDS(plot_data,file = '~/SWORD/作图/Figure2_箱线/data/plot_data_test12.rds')




# plot_data <- readRDS("~/SWORD/作图/Figure2_箱线/data/plot_data_test12.rds")
library(ggh4x)
library(ggplot2)
library(reshape2)
library(tidyverse)


for (i in 1:nrow(plot_data)) {
  if(plot_data[i,'datasets'] %in% c("deone_kidney","deone_PFC")){
    plot_data[i,'test']='test1'
  }else{
    plot_data[i,'test']='test2'
  }
}

for (i in 1:nrow(plot_data)) {
  if(plot_data[i,'datasets'] == "deone_kidney"){
    plot_data[i,'datasets2']='ER'
  }
  if(plot_data[i,'datasets'] == "deone_PFC"){
    plot_data[i,'datasets2']='PBMC'
  }
  if(plot_data[i,'datasets'] == "cross_PFC"){
    plot_data[i,'datasets2']='CRC'
  }
  if(plot_data[i,'datasets'] == "cross_PBMC"){
    plot_data[i,'datasets2']='TNBC'
  }
}

plot_data[,'methods']=factor(plot_data[,'methods'],levels=c("DWLS_sep",'DWLS_com',"BayesPrism","MuSiC","Bisque","SCDC"))
vec_color <-c('#a52a2a','#F58840','#008000','#ffa500',"#761c77",'#4861a3')


###绘图
strip = strip_nested(background_x = elem_list_rect(fill = c(rep('grey',2),
                                                            rep(c('#F58840','#9AE66E'),2)),color = rep('white',6)))

ggplot(plot_data,aes(x = methods,y = RMSE,fill = methods)) +
  geom_boxplot(aes(color = methods),#这里的fill如果不设就是空心的
               width = .6, size = .7, alpha = .5, outlier.size = 1, ) +
  theme_classic() +
  scale_color_manual(values= vec_color) +
  scale_fill_manual(values= vec_color) +
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
  coord_cartesian(ylim = c(0, 0.2))+
  # facet_wrap(~test+datasets2, ncol = 4, scales = "free")
  facet_nested(~test+datasets2,scales = 'free_y',axes = 'all',
               strip  = strip)

