library(survival)
library(survminer)
prd_after <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/GSE58812/prd_sep_init.rds")
Clinical<- readRDS("~/SWORD/应用/TNBC/GSE58812/GSE58812_metadata.rds")


# "T_cells_c3_CD4__Tfh_CXCL13"   T_cells_c8_CD8__LAG3
# "PVL_Differentiated_s3"          "PVL_Immature_s1"  "PVL_Immature_s2"

gene_exp2<-as.data.frame(t(prd_after['T_cells_c3_CD4__Tfh_CXCL13',]))
# gene_exp2<-as.data.frame(t(pvl_all))

# # #                          
# gene_exp2<-as.data.frame(t(prd_after['PVL_Immature_s1',]/pvl_all))
gene_exp<-as.data.frame(gene_exp2)

##<=的话小于的被cut掉，成为TRUE
cutoff= (gene_exp2[,1]<=quantile(as.numeric(gene_exp2[,1]))[3])#以中位数为界分为高表达及低表达组

# cutoff= (gene_exp2[,1]<=median(as.numeric(gene_exp2[,1])))#以中位数为界分为高表达及低表达组
# cutoff2=(gene_exp1[,1]>=median(as.numeric(gene_exp1[,1]))) &  (gene_exp2[,1]>=median(as.numeric(gene_exp2[,1])))#以中位数为界分为高表达及低表达组

# cutoff=gene_exp[,1]<=median(as.numeric(gene_exp[,1]))



surtime<-as.numeric(Clinical[row.names(gene_exp),'os'])
surstat<-as.numeric(Clinical[row.names(gene_exp),'death'])
# fit <- survfit(Surv(surtime, surstat) ~cutoff, data = gene_exp) # survf

fit2 <- survfit(Surv(surtime, surstat) ~  cutoff, data = gene_exp )

surv_diff <- survdiff(Surv(surtime, surstat) ~  cutoff, data = gene_exp )
# p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)#提p值

ak<-paste(cutoff)
gene_exp[,1]<-as.factor(ak)
gene_exp[,2]<-surtime
gene_exp[,3]<-surstat
colnames(gene_exp)<-c('ak2','surtime','surstat')

pairwise_survdiff(Surv(surtime, surstat) ~ ak2, data = gene_exp)
# p.value  <- pairwise_survdiff(Surv(surtime, surstat) ~ ak2, data = gene_exp)

##510 500/4 3
# GS<-
  ggsurvplot(fit2,
               pval = TRUE, pval.method = F, pval.coord= c(0.05, 0.05), pval.size = 5,#p值参数（坐标/大小等）
               legend.labs=c('High','Low'),
               conf.int = TRUE, #是否显示置信区间
               risk.table = FALSE, # 是否添加风险表
               legend = c(0.9,0.15),#c(0.8, 0.9), # legend位置坐标
               legend.title=element_blank(),  #改图例名称
               # palette="lancet", #柳叶刀配色
               # legend.title = "Gene39",legend.labs = c("High", "Low"),font.legend = 14,#legend参数调整
               # font.main = c(12, "bold", "darkblue"),font.tickslab = 12,#字体调整
               ylab= "Overall Survival",#font.x = 16, font.y =16, # 坐标轴的标题
               # ggtheme = theme_survminer(), # 更改ggplot2的主题
               palette = c("#0074b3", "#982b2b","green3"), #定义颜色,
               # legend.labs = c("normal", "high",'low'),    # 图例标签
               # risk.table = F, # 是否添加风险表
               risk.table.height = 0.3,
               # ylab=NULL
               ncensor.plot = FALSE,
               ncensor.plot.height = 0.2,
               title = 'T CD4-CXCL13',
               ggtheme = theme(
                 panel.background = element_blank(),
                 axis.text.x = element_text(family = "Times", #face = "italic",
                                            #colour = "darkred",
                                            size = rel(1.6)),
                 axis.text.y = element_text(family = "Times", #face = "italic",
                                            #colour = "darkred",
                                            size = rel(1.6)),
                 legend.text = element_text(family = "Times", #face = "italic",
                                            #colour = "darkred",
                                            size = rel(1.3)),
                 legend.title = element_text(family = "Times", #face = "italic",
                                             #colour = "darkred",
                                             size = rel(1.3)),
                 panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),
                 plot.title = element_text(family = "Times",size = 18, face = "bold",hjust=0.5),
                 axis.title = element_blank()   
               ),
               fontsize=6
)
###4  4.5