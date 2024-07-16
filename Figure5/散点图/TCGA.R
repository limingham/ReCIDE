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
colnames(data_matrix)[5]='death'
# data_matrix[,'CX_PVL']<-data_matrix[,'CXCL13']
data_matrix[,'CX_PVL']<-data_matrix[,'PVL']/data_matrix[,'PVLs3']
# data_matrix[,'CX_PVL']<-data_matrix[,'CXCL13']*data_matrix[,'PVLs1']/data_matrix[,'PVL']/data_matrix[,'PVLs3']

# /data_matrix[,'PVL']
data_matrix=data_matrix[!(data_matrix[,'death']=='0' & data_matrix[,'surtime']<730),]


options(repr.plot.width=7,repr.plot.height=7)
ggplot(data_matrix, aes(surtime, CX_PVL)) +
  geom_point(aes(fill = death, color = death, shape = death), size = 2, alpha = 0.7)+
  # geom_smooth(method = 'lm', se = FALSE, col = 'black')+
  # stat_cor(label.x = 0.02,label.y = 0.1, size=5)+
  geom_vline(xintercept=1825,linetype=3)+
  geom_hline(yintercept = 0.043,linetype=3)+
  # geom_hline(yintercept = 0.08,linetype=3)+
  scale_color_manual(values = c("#9d4649", "#348c4a"))+
  scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous( limits = c(0, 5500))+
  labs(y="PVL", x="T CD4-CXCL13")+
  ggtitle("PVL vs T CD4-CXCL13")+
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
  guides(color = guide_legend(override.aes = list(size = 4)))
###5 3.8




data_matrix1=data_matrix[data_matrix[,'CXCL13']>8e-04,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'CXCL13']<8e-04,]
table(data_matrix2[,5])


data <- matrix(c(47,23, 11, 9), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)



###############PVLs3
###############


seq(0,0.045,0.0001)[90]
# proportion= 0.0089
# p-value = 0.08431
data_matrix1=data_matrix[data_matrix[,'PVLs3']>0.0089,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'PVLs3']<0.0089,]
table(data_matrix2[,5])


data <- matrix(c(38,27, 20, 5), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)



###############model
###############



# proportion= 0.043
# p-value = 0.003135
data_matrix1=data_matrix[data_matrix[,'CX_PVL']>0.043,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'CX_PVL']<0.043,]
table(data_matrix2[,5])


data <- matrix(c(28,6, 30, 26), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)


###############PVL
###############


data_matrix1=data_matrix[data_matrix[,'PVL']>0.034,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'PVL']<0.034,]
table(data_matrix2[,5])


data <- matrix(c(44,30, 14, 2), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)




###############PVL/PVLs3
###############


data_matrix1=data_matrix[data_matrix[,'CX_PVL']>8,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'CX_PVL']<8,]
table(data_matrix2[,5])


data <- matrix(c(5,1, 53, 31), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)

