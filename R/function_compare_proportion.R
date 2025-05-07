compare_group_proportion <- function(results_list,dir_diff_celltype_output){
  df_category1 = as.data.frame(t(results_list[[1]]))
  df_category2 = as.data.frame(t(results_list[[2]]))
  
  
  col_inter = intersect(colnames(df_category1),colnames(df_category2))
  diff1 <- setdiff(colnames(df_category1), col_inter)
  diff2 <- setdiff(colnames(df_category2), col_inter)
  if(length(diff1) > 0){df_category2[,diff1] = 0}
  
  if(length(diff2) > 0){df_category1[,diff2] = 0}
  
  df_category1 = df_category1[,sort(colnames(df_category1))]
  df_category2 = df_category2[,sort(colnames(df_category2))]
  
  p_values <- numeric(length = ncol(df_category1))
  trends <- character(length = ncol(df_category2))
  diff <- numeric(length = ncol(df_category2))
  
  for (j in seq_along(colnames(df_category1))) {
    col_name <- colnames(df_category1)[j]
    
    df_category1_in = df_category1[, col_name]
    
    
    df_category2_in = df_category2[, col_name]
    
    
    pvalue <- wilcox.test(df_category1_in, df_category2_in, alternative = "two.sided")$p.value
    p_values[j] <- pvalue
    diff[j] <- abs(median(df_category1_in)-median(df_category2_in))
    if (diff[j]==0){
      diff[j] <- abs(mean(df_category1_in)-mean(df_category2_in))
    }
    #   
    if (median(df_category1_in)<median(df_category2_in)) {
      trends[j] <- '+'  
    } else if (median(df_category1_in)>median(df_category2_in)) {
      trends[j] <- '-'
    } else {
      trends[j] <- '0'  
    }
  }
  adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
  
  re_df <- data.frame(
    Celltype = colnames(df_category1),
    PValue = adjusted_p_values,
    Trend = trends,
    Diff=diff
  )
  saveRDS(re_df,file=dir_diff_celltype_output)
}






