library(ggplot2)
library(ggpubr)

prd_after <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/所有显著细胞类型箱线图/GSE164458/prd_sep_40.rds")
prd_after=prd_after[sort(row.names(prd_after)),]

prd_after=as.data.frame(t(prd_after))

# prd_after=subset(prd_after,T_cells_c3_CD4__Tfh_CXCL13 >5.192176e-03 | T_cells_c3_CD4__Tfh_CXCL13 <1.452157e-03)

pdata<- readRDS("~/SWORD/应用/TNBC/GSE164458/GSE164458_metadata.rds")
row.names(pdata)<-pdata[,2]

# pdata=pdata[row.names(prd_after),]

all(row.names(prd_after)==row.names(pdata))
# pdata<-pdata[row.names(prd_after),]
# pdata<-pdata[row.names(prd_after),]

# "T_cells_c3_CD4__Tfh_CXCL13"
# "PVL_Differentiated_s3"          "PVL_Immature_s1"  "PVL_Immature_s2"


p_list=list()

for (k in 1:ncol(prd_after)) {


boxdata<-as.data.frame(prd_after[,k])


row.names(boxdata)=row.names(prd_after)

# boxdata<-as.data.frame(t(prd_after['B_cell',]))
boxdata[,2]<-pdata[,'pathologic_complete_response']


boxdata[,1]<-as.numeric(boxdata[,1])

colnames(boxdata)<-c("X1","X2")

library(dplyr)


y_high1=quantile(boxdata[,1],probs =0.99)
y_high2=quantile(boxdata[,1],probs =0.97)

p_list[[k]]=ggplot(boxdata, aes(x = X2, y = X1)) +
  geom_boxplot(aes(color = X2),#这里的fill如果不设就是空心的
               width = .6, size = .7, alpha = .5, outlier.size = .1, ) +
  stat_boxplot(geom = "errorbar",
               aes(ymin = ..ymax.., color = X2),
               #这是调上界长度和粗细
               width = 0.5, size = 1) +
  stat_boxplot(geom = "errorbar",
               #这是调下界长度和粗细
               aes(ymax = ..ymin.., color = X2), 
               width = 0.5, size = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, face = "bold", angle = 50, hjust = 1),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
  )+
  # labs (x = "", y = "") +
  #labs (x = "Type", y = "Accuracy Rate(%)", title = "Total contacts", subtitle = "(chr3)") +
  #scale_fill_manual(values = c("grey91","thistle1","blanchedalmond","darkseagreen1","azure2")) +
  scale_fill_manual(values = c("#982b2b","#0074b3","red3")) +
  # scale_color_manual(values= c("cadetblue1","powderblue","palegreen4","peachpuff1","pink2","red3")) +
  # scale_fill_manual(values = c("powderblue","pink2")) +
  scale_color_manual(values = c("#982b2b","#0074b3","red3")) +
  guides(fill = "none", color = "none")+
  geom_point(aes(colour = factor(X2)),alpha=0.4)+
  theme(panel.background =element_blank(),
        axis.text.x = element_text(family = "Times", #face = "italic",
                                   #colour = "darkred",
                                   size = rel(1.4)),
        axis.text.y = element_text(family = "Times", #face = "italic",
                                   #colour = "darkred",
                                   size = rel(1.4)),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        plot.title = element_text(family = "Times",size = 12),
        axis.title = element_blank()
  )+  
  coord_cartesian(ylim = c(0, y_high1))+
  stat_compare_means(label.y=y_high2)+
  ggtitle(colnames(prd_after)[k])
  ###3.2 3

}
names(p_list)=colnames(prd_after)



library(cowplot)
library(gridExtra)
# p_list[[1]]+p_list[[2]]+p_list[[3]]+p_list[[4]]+p_list[[5]]+p_list[[6]]+p_list[[7]]+
# p_list[[8]]+p_list[[9]]+p_list[[10]]+p_list[[11]]+p_list[[12]]+p_list[[13]]+p_list[[14]]+
# p_list[[15]]+p_list[[16]]+p_list[[17]]+p_list[[18]]+p_list[[19]]+p_list[[20]]+p_list[[21]]+
# p_list[[22]]+p_list[[23]]+p_list[[24]]+p_list[[25]]+p_list[[26]]+p_list[[27]]+p_list[[28]]+
# p_list[[29]]+p_list[[30]]+p_list[[31]]+p_list[[32]]+p_list[[33]]+p_list[[34]]+p_list[[35]]+
# p_list[[36]]+p_list[[37]]+p_list[[38]]+p_list[[39]]+p_list[[40]]+p_list[[41]]+p_list[[42]]+
# p_list[[43]]+p_list[[44]]+p_list[[45]]+p_list[[46]]+p_list[[47]]+p_list[[48]]+p_list[[49]]
#   
# plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]])
# grid.arrange(p_list, print = TRUE,  
#                     ncol = 4, nrow = 13) #定义行数和列数
# 
# 
grid.arrange(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],p_list[[5]],p_list[[6]],p_list[[7]],
          p_list[[8]],p_list[[9]],p_list[[10]],p_list[[11]],p_list[[12]],p_list[[13]],p_list[[14]],
          p_list[[15]],p_list[[16]],p_list[[17]],p_list[[18]],p_list[[19]],p_list[[20]],p_list[[21]],
          p_list[[23]],p_list[[24]],p_list[[25]],p_list[[26]],p_list[[28]],
          p_list[[29]],p_list[[30]],p_list[[32]],p_list[[34]],p_list[[35]],
          p_list[[36]],p_list[[37]],p_list[[38]],p_list[[40]],p_list[[41]],
          p_list[[43]],p_list[[44]],p_list[[46]],p_list[[47]],p_list[[48]],p_list[[49]],nrow=7,ncol=6)
###4  4.5
# ####14,11


# grid.arrange(p_list[[1]],p_list[[13]],p_list[[25]],p_list[[37]],
#              p_list[[2]],p_list[[14]],p_list[[26]],p_list[[38]],
#              p_list[[3]],p_list[[15]],
#              p_list[[4]],p_list[[16]], p_list[[28]],p_list[[40]],
#              p_list[[5]],p_list[[17]],p_list[[29]],p_list[[41]],
#              p_list[[6]],p_list[[18]],p_list[[30]],
#              p_list[[7]],p_list[[19]],p_list[[43]],
#              p_list[[8]],p_list[[20]],p_list[[32]],p_list[[44]],
#              p_list[[9]],p_list[[21]],
#              p_list[[10]],p_list[[34]],p_list[[46]],
#              p_list[[11]],p_list[[23]],p_list[[35]],p_list[[47]],
#              p_list[[12]],p_list[[24]],p_list[[36]],p_list[[48]],
#              p_list[[49]],nrow=7,ncol=6)
###10  40