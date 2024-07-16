library(ggplot2)
library(ggpubr)

prd_after <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/GSE164458/prd_sep_40.rds")
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
pvl_all=prd_after[,'PVL_Immature_s1']/(prd_after[,'PVL_Differentiated_s3']+prd_after[,'PVL_Immature_s1']+prd_after[,'PVL_Immature_s2']+prd_after[,'Cycling_PVL'])
 # prd_after[,'PVL_Differentiated_s3']+prd_after[,'PVL_Immature_s1']+prd_after[,'PVL_Immature_s2']+prd_after[,'Cycling_PVL']
# boxdata<-as.data.frame(prd_after[,'PVL_Immature_s1']/pvl_all)
boxdata<-as.data.frame(pvl_all)


# boxdata<-as.data.frame(prd_after[,'T_cells_c3_CD4__Tfh_CXCL13'])


row.names(boxdata)=row.names(prd_after)

# boxdata<-as.data.frame(t(prd_after['B_cell',]))
boxdata[,2]<-pdata[,'pathologic_complete_response']

# row.names(boxdata)<-row.names(prd_after)
# boxdata<-subset(boxdata,boxdata[,2] %in% c('status: CTL','status: AD'))
# boxdata<-subset(boxdata,pdata[row.names(boxdata),"included in case -control study:ch1"] %in% 'yes')


# boxdata<-subset(boxdata,boxdata[,2] %in% c('status: MCI','status: CTL'))
boxdata[,1]<-as.numeric(boxdata[,1])

colnames(boxdata)<-c("X1","X2")

library(dplyr)

# boxdata<-P_sep[["data"]]
# boxdata=boxdata %>% mutate(X2=case_when(X2 == "status: CTL" ~ "CTL",
#                                         X2 == "status: MCI" ~ "MCI",
#                                         X2 == "status: AD" ~ "AD"))

# P_sep<-
ggplot(boxdata, aes(x = X2, y = X1)) +
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
  # coord_cartesian(ylim = c(0, 0.12))+
  stat_compare_means(label.y=0.75)
  ###3.5 3
