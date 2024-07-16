library(survival)
library(survminer)
prd_after <- readRDS("~/ReCIDE/应用_TNBC/论文中结果/GSE58812/prd_sep_init.rds")
Clinical<- readRDS("~/ReCIDE/应用_TNBC/GSE58812/GSE58812_metadata.rds")

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
  scale_y_continuous(limits = c(0, 1.3))+
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

###############CXCL13
###############
FT=c()
k=1
for (i in seq(0,2.2e-02,0.0001)) {
  data_matrix1=data_matrix[data_matrix[,'CXCL13']>i,]
  data1=table(data_matrix1[,5])
  
  data_matrix2=data_matrix[data_matrix[,'CXCL13']<i,]
  data2=table(data_matrix2[,5])
  
  
  data <- matrix(c(data1[1], data1[2], data2[1], data2[2]), nrow = 2) 
  colnames(data) <- c("组1", "组2") 
  rownames(data) <- c("暴露", "未暴露")
  
  data[is.na(data)]=0
  # 执行费舍尔精确检验 
  FT[k] =  fisher.test(data)  
  k=k+1
}


order(as.numeric(FT))
seq(0,2.2e-02,0.0001)[9]
# proportion= 8e-04
# p-value = 0.00177
data_matrix1=data_matrix[data_matrix[,'CXCL13']>8e-04,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'CXCL13']<8e-04,]
table(data_matrix2[,5])


data <- matrix(c(69,19, 6, 10), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)



###############PVLs3
###############

FT=c()
k=1
for (i in seq(0,0.045,0.0001)) {
  data_matrix1=data_matrix[data_matrix[,'PVLs3']>i,]
  data1=table(data_matrix1[,5])
  
  data_matrix2=data_matrix[data_matrix[,'PVLs3']<i,]
  data2=table(data_matrix2[,5])
  
  
  data <- matrix(c(data1[1], data1[2], data2[1], data2[2]), nrow = 2) 
  colnames(data) <- c("组1", "组2") 
  rownames(data) <- c("暴露", "未暴露")
  
  data[is.na(data)]=0
  # 执行费舍尔精确检验 
  FT[k] =  fisher.test(data)  
  k=k+1
}


min(as.numeric(FT))
order(as.numeric(FT))

seq(0,0.045,0.0001)[90]
# proportion= 0.0089
# p-value = 0.004548
data_matrix1=data_matrix[data_matrix[,'PVLs3']>0.0089,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'PVLs3']<0.0089,]
table(data_matrix2[,5])


data <- matrix(c(28,20, 47, 9), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)



###############model
###############

FT=c()
k=1
for (i in seq(0,0.2,0.001)) {
  data_matrix1=data_matrix[data_matrix[,'CX_PVL']>i,]
  data1=table(data_matrix1[,5])
  
  data_matrix2=data_matrix[data_matrix[,'CX_PVL']<i,]
  data2=table(data_matrix2[,5])
  
  
  data <- matrix(c(data1[1], data1[2], data2[1], data2[2]), nrow = 2) 
  colnames(data) <- c("组1", "组2") 
  rownames(data) <- c("暴露", "未暴露")
  
  data[is.na(data)]=0
  # 执行费舍尔精确检验 
  FT[k] =  fisher.test(data)  
  k=k+1
}


min(as.numeric(FT))
order(as.numeric(FT))

seq(0,0.2,0.001)[44]
# proportion= 0.043
# p-value = 0.003135
data_matrix1=data_matrix[data_matrix[,'CX_PVL']>0.043,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'CX_PVL']<0.043,]
table(data_matrix2[,5])


data <- matrix(c(55,12, 20, 17), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)





###############PVLs1/PVL
###############

data_p=data_matrix[,'PVLs1']/data_matrix[,'PVL']

FT=c()
k=1
for (i in seq(0,0.7,0.01)) {
  data_matrix1=data_matrix[data_p>i,]
  data1=table(data_matrix1[,5])
  
  data_matrix2=data_matrix[data_p<i,]
  data2=table(data_matrix2[,5])
  
  
  data <- matrix(c(data1[1], data1[2], data2[1], data2[2]), nrow = 2) 
  colnames(data) <- c("组1", "组2") 
  rownames(data) <- c("暴露", "未暴露")
  
  data[is.na(data)]=0
  # 执行费舍尔精确检验 
  FT[k] =  fisher.test(data)  
  k=k+1
}


min(as.numeric(FT))
order(as.numeric(FT))

seq(0,0.7,0.01)[39]
# proportion= 0.043
# p-value = 0.003135
data_matrix1=data_matrix[data_p>0.38,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_p<0.38,]
table(data_matrix2[,5])


data <- matrix(c(31,6, 44, 23), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)



##############PVL
###############
FT=c()
k=1
for (i in seq(0,0.1,0.001)) {
  data_matrix1=data_matrix[data_matrix[,'PVL']>i,]
  data1=table(data_matrix1[,5])
  
  data_matrix2=data_matrix[data_matrix[,'PVL']<i,]
  data2=table(data_matrix2[,5])
  
  
  data <- matrix(c(data1[1], data1[2], data2[1], data2[2]), nrow = 2) 
  colnames(data) <- c("组1", "组2") 
  rownames(data) <- c("暴露", "未暴露")
  
  data[is.na(data)]=0
  # 执行费舍尔精确检验 
  FT[k] =  fisher.test(data)  
  k=k+1
}


min(as.numeric(FT))
order(as.numeric(FT))
seq(0,0.1,0.001)[35]
# proportion= 0.0337
# p-value = 0.00177
data_matrix1=data_matrix[data_matrix[,'PVL']>0.034,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'PVL']<0.034,]
table(data_matrix2[,5])


data <- matrix(c(31,19, 44, 10), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)




##############PVLs3/PVL
###############
FT=c()
k=1
for (i in seq(1,18,0.1)) {
  data_matrix1=data_matrix[data_matrix[,'CX_PVL']>i,]
  data1=table(data_matrix1[,5])
  
  data_matrix2=data_matrix[data_matrix[,'CX_PVL']<i,]
  data2=table(data_matrix2[,5])
  
  
  data <- matrix(c(data1[1], data1[2], data2[1], data2[2]), nrow = 2) 
  colnames(data) <- c("组1", "组2") 
  rownames(data) <- c("暴露", "未暴露")
  
  data[is.na(data)]=0
  # 执行费舍尔精确检验 
  FT[k] =  fisher.test(data)  
  k=k+1
}


min(as.numeric(FT))
order(as.numeric(FT))

seq(1,18,0.1)[71]
# proportion= 0.043
# p-value = 0.003135
data_matrix1=data_matrix[data_matrix[,'CX_PVL']>8,]
table(data_matrix1[,5])

data_matrix2=data_matrix[data_matrix[,'CX_PVL']<8,]
table(data_matrix2[,5])


data <- matrix(c(16,1, 59, 28), nrow = 2) 
colnames(data) <- c("组1", "组2") 
rownames(data) <- c("暴露", "未暴露")

# 执行费舍尔精确检验 
# chisq.test(data)
fisher.test(data)
