library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
# prd_bisque_df <- readRDS("~/ReCIDE/benchmark测试/deone_kidney/bisque/prd_bisque_df.rds")
# prd_bisque_df['NEU',]=0
# saveRDS(prd_bisque_df,file="~/ReCIDE/benchmark测试/deone_kidney/bisque/prd_bisque_df.rds")
key_df <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/EXP_and_KEY/KEY.rds")
row.names(key_df) <- str_replace(row.names(key_df), "\\+", "_")
row.names(key_df) <- str_replace(row.names(key_df), "-", "_")
row.names(key_df) <- str_replace(row.names(key_df), " ", "_")
row.names(key_df) <- str_replace(row.names(key_df), " ", "_")
row.names(key_df) <- str_replace(row.names(key_df), " ", "_")

key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]

prd_DWLS_com <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/DWLS_com/prd_df.rds")
prd_DWLS_com[is.na(prd_DWLS_com)]=0

# prd_DWLS_com_more=row.names(prd_DWLS_com)[!(row.names(prd_DWLS_com) %in% row.names(key_df))]
# key_df_more=row.names(key_df)[!(row.names(key_df) %in% row.names(prd_DWLS_com))]

# prd_DWLS_com[key_df_more,]=0
# key_df[prd_DWLS_com_more,]=0


# prd_BayesPrism <- readRDS("~/ReCIDE/benchmark测试/high_res_TNBC/bayesprism/prd_df.rds")
prd_BayesPrism <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bayesprism/prd_df.rds")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), "\\+", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), "-", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), " ", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), " ", "_")
row.names(prd_BayesPrism) <- str_replace(row.names(prd_BayesPrism), " ", "_")
prd_BayesPrism[is.na(prd_BayesPrism)]=0

# prd_BayesPrism[is.na(prd_BayesPrism)]=0
# prd_BayesPrism[prd_DWLS_com_more,]=0


prd_Bisque<- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/bisque/prd_df.rds")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), "\\+", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), "-", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), " ", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), " ", "_")
row.names(prd_Bisque) <- str_replace(row.names(prd_Bisque), " ", "_")
prd_Bisque[is.na(prd_Bisque)]=0

# prd_Bisque[is.na(prd_Bisque)]=0
# prd_Bisque[prd_DWLS_com_more,]=0



prd_DWLS_sep <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/DWLS_sep/prd_df.rds")
# prd_DWLS_sep[prd_DWLS_com_more,]=0

prd_MuSiC  <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/music/prd_df.rds")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), "\\+", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), "-", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), " ", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), " ", "_")
row.names(prd_MuSiC) <- str_replace(row.names(prd_MuSiC), " ", "_")
prd_MuSiC[is.na(prd_MuSiC)]=0

key_df_more=row.names(key_df)[!(row.names(key_df) %in% row.names(prd_MuSiC))]
# prd_MuSiC[is.na(prd_MuSiC)]=0
prd_MuSiC[key_df_more,]=0



prd_SCDC <- readRDS("~/ReCIDE/benchmark测试_最终/high_res_TNBC/scdc/prd_df.rds")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), "\\+", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), "-", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), " ", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), " ", "_")
row.names(prd_SCDC) <- str_replace(row.names(prd_SCDC), " ", "_")
prd_SCDC[is.na(prd_SCDC)]=0

# prd_SCDC[is.na(prd_SCDC)]=0
# prd_SCDC[prd_DWLS_com_more,]=0


prd_list=list(prd_BayesPrism,
              prd_Bisque,
              prd_DWLS_com,
              prd_DWLS_sep,
              prd_MuSiC,
              prd_SCDC)

names(prd_list)=c('BayesPrism',
                  'Bisque',
                  'DWLS_com',
                  'DWLS_sep',
                  'MuSiC',
                  'SCDC')

for (i in 1:length(prd_list)) {
  prd_list[[i]]=prd_list[[i]][sort(row.names(prd_list[[i]])),sort(colnames(prd_list[[i]]))]
}



key_df=key_df[sort(row.names(key_df)),sort(colnames(key_df))]

for (k in 1:length(prd_list)) {
  print(row.names(prd_list[[k]])==row.names(key_df))
}


data_mean = apply(key_df,1,mean)

names_mean=names(data_mean[data_mean<0.02])

key_df=key_df[names_mean,]

# key_df2=apply(key_df,1,scale)
# 
# key_df2=as.data.frame(t(key_df2))
# colnames(key_df2)=colnames(key_df)


for (k in 1:length(prd_list)) {
  prd_list[[k]]=prd_list[[k]][names_mean,]
}

# RMSE_list=list()
# for (i in 1:length(prd_list)) {
#   df_in=c()
#   for (j in 1:ncol(prd_list[[i]])) {
#     df_in=c(df_in,ModelMetrics::rmse(prd_list[[i]][,j],key_df[,j]))
#   }#ModelMetrics::rmse
#   df_in=as.data.frame(df_in)
#   df_in[,2]=names(prd_list)[i]
#   colnames(df_in)=c('RMSE','method')
#   RMSE_list[[i]]=df_in
#   names(RMSE_list)[i]=names(prd_list)[i]
# }
RMSE_list=list()
for (i in 1:length(prd_list)) {
  df_in=c()
  for (j in 1:ncol(prd_list[[i]])) {
    # prd_list[[i]][,j]= prd_list[[i]][,j]/sum(prd_list[[i]][,j])*sum(key_df[,j])
    # prd_list[[i]][,j]= prd_list[[i]][,j]
    
    df_in=c(df_in,ModelMetrics::rmse(as.numeric(prd_list[[i]][,j]),as.numeric(key_df[,j])))
    # df_in=c(df_in,ModelMetrics::rmse(as.numeric(prd_list[[i]][,j]),as.numeric(key_df[,j])))
    
    }#ModelMetrics::rmse
  df_in=as.data.frame(df_in)
  df_in[,2]=names(prd_list)[i]
  colnames(df_in)=c('RMSE','method')
  RMSE_list[[i]]=df_in
  names(RMSE_list)[i]=names(prd_list)[i]
}



plot_data=RMSE_list[[1]]

for (i in 2:length(RMSE_list)) {
  plot_data=rbind(plot_data,
                  RMSE_list[[i]])
}

df_median=as.data.frame(group_by(plot_data,method) %>% summarize_each(median))

df_median<-arrange(df_median, df_median[,'RMSE'])




plot_data[,'method']=factor(plot_data[,'method'],levels=c("DWLS_sep","SCDC","DWLS_com",
                                                            'BayesPrism',"MuSiC","Bisque"))


vec_color <-c('#a52a2a','#4861a3','#F58840','#008000','#ffa500',"#761c77")
library(RColorBrewer)



ggplot(plot_data, aes(x= reorder(method,RMSE), y=RMSE,fill=method)) + 
  geom_boxplot(aes(color = method),#这里的fill如果不设就是空心的
               size = .7, alpha = .5, outlier.size = 1,
               position = position_dodge(width = 0.5), width = 0.5) +
  # position_dodge是箱子间距离  width是箱宽
  # geom_boxplot(data=subset(plot_data, methods %in% 'empty'), width = .5) +
  theme_classic() +
  geom_abline(intercept = median(df_median[1,2]), slope = 0, col = 'Black',linewidth=0.45,linetype = 2)+
  theme(
    axis.text.x = element_text(size = 10, face = "plain", angle = -45),
    axis.text.y = element_text(size = 10, face = "plain"),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title = element_text(size = 10, face = "plain"),
    plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "plain", hjust = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    # legend.position = 'none'
    # axis.title = element_text(size = 8)
  )+
  labs (x = "", y = "") +
  scale_color_manual(values= vec_color) +
  # guides(fill = "none", color = "none")+
  scale_fill_manual(values= vec_color) +
  coord_cartesian(ylim = c(0, 0.07))+
  # coord_cartesian(ylim = c(0.01, 0.13))+
  labs (x = "", y = "")#+

##4.5 3