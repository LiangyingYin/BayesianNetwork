#add some rank inforamtion
Fisher_IDA_impExpr <- function(dt_IDA, topK_percent= c(0.1,0.2,0.3,0.4,0.5) ) {
  
  res = matrix(nrow= length(topK_percent), ncol=3) 
  
  for (i in 1:length(topK_percent) ) {
    
    topK_num = round(topK_percent[i]*genes_num) 
    ncol_IDA = ncol(dt_IDA)
    rank_true_IDA = rank(-abs(dt_IDA$true_IDA))
    rank_estimate_IDA = rank(-abs(dt_IDA$Gene_imputedata_regression_Cof))

    
    topK_true = as.numeric(rank_true_IDA<= topK_num) 
    topK_est = as.numeric(rank_estimate_IDA<= topK_num) 
    
    final_table = table(topK_true, topK_est)	
    final_table
    
    if ((nrow(final_table)>=2) && (ncol(final_table)>=2))
    {
      fisher_test_obj = fisher.test(final_table, alternative="greater")
      res[i, ] = c(topK_percent[i], fisher_test_obj$estimate, fisher_test_obj$p.value)	
    }
    
  }
  colnames(res) <- c("topK_percent", "OR", "pval_IDA")
  return(res)
}   

Fisher_IDA <- function(dt_IDA, topK_percent= c(0.1,0.2,0.3,0.4,0.5) ) {
  
  res = matrix(nrow= length(topK_percent), ncol=3) 
  
  for (i in 1:length(topK_percent) ) {
    
    topK_num = round(topK_percent[i]*genes_num)
    ncol_IDA = ncol(dt_IDA)
    rank_true_IDA = rank(-abs(dt_IDA$true_IDA))
    rank_estimate_IDA = rank(-abs(dt_IDA$Estimate_IDA_final))
    
    
   
    
    topK_true = as.numeric(rank_true_IDA<= topK_num) 
    topK_est = as.numeric( rank_estimate_IDA<= topK_num) 
    
    final_table = table(topK_true, topK_est)	
    fisher_test_obj = fisher.test(final_table, alternative="greater")
    res[i, ] = c(topK_percent[i], fisher_test_obj$estimate, fisher_test_obj$p.value)		
  }
  colnames(res) <- c("topK_percent", "OR", "pval_IDA") 
  return(res)
}    

Fisher_jointIDA_impExpr <- function(dt_IDA, topK_percent= c(0.1,0.2,0.3,0.4,0.5) ) {
  res_joint = matrix(nrow= length(topK_percent), ncol=3)   
  
  for (i in 1:length(topK_percent) ) {
    
    topK_num = round(topK_percent[i]*genes_num)
    ncol_IDA = ncol(dt_IDA)
    rank_true_IDA = rank(-abs(dt_IDA$true_jointIDA))
    rank_estimate_IDA = rank(-abs(dt_IDA$Gene_imputedata_regression_Cof ))
    
    topK_true = as.numeric(rank_true_IDA<= topK_num) 
    topK_est = as.numeric( rank_estimate_IDA<= topK_num) 
    
    final_table = table(topK_true, topK_est)
    fisher_test_obj = fisher.test(final_table, alternative="greater")
    
    res_joint[i, ] = c(topK_percent[i], fisher_test_obj$estimate, fisher_test_obj$p.value)
    
  }
  colnames(res_joint) <- c("topK_percent", "OR", "pval_JointIDA") 
  return(res_joint)
}    

Fisher_jointIDA <- function(dt_IDA, topK_percent= c(0.1,0.2,0.3,0.4,0.5) ) {
  res_joint = matrix(nrow= length(topK_percent), ncol=3)   
  
  for (i in 1:length(topK_percent) ) {
    
    topK_num = round(topK_percent[i]*genes_num)
    ncol_IDA = ncol(dt_IDA)
    rank_true_IDA = rank(-abs(dt_IDA$true_jointIDA))
    rank_estimate_IDA = rank(-abs(dt_IDA$est_jointIDA ))
    
    topK_true = as.numeric(rank_true_IDA<= topK_num) 
    topK_est = as.numeric( rank_estimate_IDA<= topK_num) 
    
    final_table = table(topK_true, topK_est)
    fisher_test_obj = fisher.test(final_table, alternative="greater")
    
    res_joint[i, ] = c(topK_percent[i], fisher_test_obj$estimate, fisher_test_obj$p.value)
    
  }
  colnames(res_joint) <- c("topK_percent", "OR", "pval_JointIDA") 
  return(res_joint)
}    




# updated by fyn 9/2/2021. Comparsion between imputedcov and rawcov
Fisher_IDA_Forcomparsion <- function(dt_IDA, topK_percent= c(0.1,0.2,0.3,0.4,0.5) ) {
  
  res = matrix(nrow= length(topK_percent), ncol=3) 
  
  for (i in 1:length(topK_percent) ) {
    
    topK_num = round(topK_percent[i]*genes_num)
    ncol_IDA = ncol(dt_IDA)
    rank_true_IDA = rank(-abs(dt_IDA$true_IDA))
    rank_estimate_IDA = rank(-abs(dt_IDA$Estimate_IDA))
    
  
    topK_true = as.numeric(rank_true_IDA<= topK_num) 
    topK_est = as.numeric( rank_estimate_IDA<= topK_num) 
    
    final_table = table(topK_true, topK_est)	
    fisher_test_obj = fisher.test(final_table, alternative="greater")
    res[i, ] = c(topK_percent[i], fisher_test_obj$estimate, fisher_test_obj$p.value)		
  }
  colnames(res) <- c("topK_percent", "OR", "pval_IDA") 
  return(res)
}    


# updated by fyn 9/2/2021. Comparsion between imputedcov and rawcov
Fisher_jointIDA_Forcomparsion <- function(dt_IDA, topK_percent= c(0.1,0.2,0.3,0.4,0.5) ) {
  res_joint = matrix(nrow= length(topK_percent), ncol=3)   
  
  for (i in 1:length(topK_percent) ) {
    
    topK_num = round(topK_percent[i]*genes_num)
    ncol_IDA = ncol(dt_IDA)
    rank_true_IDA = rank(-abs(as.numeric(dt_IDA$true_jointIDA)))
    rank_estimate_IDA = rank(-abs(as.numeric(dt_IDA$est_jointIDA )))
    
    topK_true = as.numeric(rank_true_IDA<= topK_num) 
    topK_est = as.numeric( rank_estimate_IDA<= topK_num) 
    
    final_table = table(topK_true, topK_est)
    fisher_test_obj = fisher.test(final_table, alternative="greater")
    
    res_joint[i, ] = c(topK_percent[i], fisher_test_obj$estimate, fisher_test_obj$p.value)
    
  }
  colnames(res_joint) <- c("topK_percent", "OR", "pval_JointIDA") 
  return(res_joint)
}    


#add by fyn 13/2/2020 caculate the average of each element in input
Fisher_Results_Average <- function(list_fisher_result)

{
  
  # list fisher_IDA_impExpr
  list_0.1percent_OR = list()
  list_0.2percent_OR = list()
  list_0.3percent_OR = list()
  list_0.4percent_OR = list()
  list_0.5percent_OR = list()
  
  list_0.1percent_pvalue = list()
  list_0.2percent_pvalue = list()
  list_0.3percent_pvalue = list()
  list_0.4percent_pvalue = list()
  list_0.5percent_pvalue = list()
  
  
  # list fisher_IDA
  list_0.1percent_fisher_IDA_OR = list()
  list_0.2percent_fisher_IDA_OR = list()
  list_0.3percent_fisher_IDA_OR = list()
  list_0.4percent_fisher_IDA_OR = list()
  list_0.5percent_fisher_IDA_OR = list()
  
  list_0.1percent_fisher_IDA_pvalue = list()
  list_0.2percent_fisher_IDA_pvalue = list()
  list_0.3percent_fisher_IDA_pvalue = list()
  list_0.4percent_fisher_IDA_pvalue = list()
  list_0.5percent_fisher_IDA_pvalue = list()
  
  
  
  # list fisher_jointIDA_impExpr
  list_0.1percent_fisher_jointIDA_impExpr_OR = list()
  list_0.2percent_fisher_jointIDA_impExpr_OR = list()
  list_0.3percent_fisher_jointIDA_impExpr_OR = list()
  list_0.4percent_fisher_jointIDA_impExpr_OR = list()
  list_0.5percent_fisher_jointIDA_impExpr_OR = list()
  
  list_0.1percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.2percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.3percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.4percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.5percent_fisher_jointIDA_impExpr_pvalue = list()
  
  
  # list fisher_jointIDA
  list_0.1percent_fisher_jointIDA_OR = list()
  list_0.2percent_fisher_jointIDA_OR = list()
  list_0.3percent_fisher_jointIDA_OR = list()
  list_0.4percent_fisher_jointIDA_OR = list()
  list_0.5percent_fisher_jointIDA_OR = list()
  
  list_0.1percent_fisher_jointIDA_pvalue = list()
  list_0.2percent_fisher_jointIDA_pvalue = list()
  list_0.3percent_fisher_jointIDA_pvalue = list()
  list_0.4percent_fisher_jointIDA_pvalue = list()
  list_0.5percent_fisher_jointIDA_pvalue = list()
  
  
  for(i in 1:length(list_fisher_result))
  {
    fisher_result = list_fisher_result[[i]]
    
    # fisher_IDA_impExpr
    fisher_IDA_impExpr = fisher_result[2:6,]
    list_0.1percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[1,2])
    list_0.2percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[2,2])
    list_0.3percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[3,2])
    list_0.4percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[4,2])
    list_0.5percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[5,2])
    list_0.1percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[1,3])
    list_0.2percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[2,3])
    list_0.3percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[3,3])
    list_0.4percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[4,3])
    list_0.5percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[5,3])
    
    
    # fisher_IDA
    fisher_IDA = fisher_result[8:12,]
    list_0.1percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[1,2])
    list_0.2percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[2,2])
    list_0.3percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[3,2])
    list_0.4percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[4,2])
    list_0.5percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[5,2])
    list_0.1percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[1,3])
    list_0.2percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[2,3])
    list_0.3percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[3,3])
    list_0.4percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[4,3])
    list_0.5percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[5,3])
    
    
    # fisher_jointIDA_impExpr
    fisher_jointIDA_impExpr = fisher_result[14:18,]
    list_0.1percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[1,2])
    list_0.2percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[2,2])
    list_0.3percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[3,2])
    list_0.4percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[4,2])
    list_0.5percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[5,2])
    list_0.1percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[1,3])
    list_0.2percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[2,3])
    list_0.3percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[3,3])
    list_0.4percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[4,3])
    list_0.5percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[5,3])
    
    
    # fisher_jointIDA
    fisher_jointIDA = fisher_result[20:24,]
    list_0.1percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[1,2])
    list_0.2percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[2,2])
    list_0.3percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[3,2])
    list_0.4percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[4,2])
    list_0.5percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[5,2])
    list_0.1percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[1,3])
    list_0.2percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[2,3])
    list_0.3percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[3,3])
    list_0.4percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[4,3])
    list_0.5percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[5,3])

    
  }
  
  
  
  fisher_IDA_impExpr_OR_average = c(mean(unlist(list_0.1percent_OR)),mean(unlist(list_0.2percent_OR)),mean(unlist(  list_0.3percent_OR)),mean(unlist(  list_0.4percent_OR)),mean(unlist(  list_0.5percent_OR )))
  fisher_IDA_impExpr_pvalue_average = c(mean(unlist(list_0.1percent_pvalue)),mean(unlist(list_0.2percent_pvalue)),mean(unlist(  list_0.3percent_pvalue)),mean(unlist(  list_0.4percent_pvalue)),mean(unlist(  list_0.5percent_pvalue)))
  
  fisher_IDA_impExpr_average = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_IDA_impExpr_OR_average,fisher_IDA_impExpr_pvalue_average)
  
  
  fisher_IDA_OR_average = c(mean(unlist(list_0.1percent_fisher_IDA_OR)), mean(unlist(  list_0.2percent_fisher_IDA_OR )), mean(unlist(  list_0.3percent_fisher_IDA_OR)), mean(unlist(  list_0.4percent_fisher_IDA_OR)), mean(unlist(  list_0.5percent_fisher_IDA_OR)))
  fisher_IDA_pvalue_average = c(mean(unlist(  list_0.1percent_fisher_IDA_pvalue)), mean(unlist(  list_0.2percent_fisher_IDA_pvalue)), mean(unlist(  list_0.3percent_fisher_IDA_pvalue)), mean(unlist(  list_0.4percent_fisher_IDA_pvalue)), mean(unlist(  list_0.5percent_fisher_IDA_pvalue)))
  fisher_IDA_average = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_IDA_OR_average,fisher_IDA_pvalue_average)
  
  fisher_jointIDA_impExpr_OR_average = c(mean(unlist(list_0.1percent_fisher_jointIDA_impExpr_OR)), mean(unlist(  list_0.2percent_fisher_jointIDA_impExpr_OR )), mean(unlist(  list_0.3percent_fisher_jointIDA_impExpr_OR)), mean(unlist(  list_0.4percent_fisher_jointIDA_impExpr_OR)), mean(unlist(  list_0.5percent_fisher_jointIDA_impExpr_OR)))
  fisher_jointIDA_impExpr_pvalue_average = c(mean(unlist(  list_0.1percent_fisher_jointIDA_impExpr_pvalue)), mean(unlist(  list_0.2percent_fisher_jointIDA_impExpr_pvalue)), mean(unlist(  list_0.3percent_fisher_jointIDA_impExpr_pvalue)), mean(unlist(  list_0.4percent_fisher_jointIDA_impExpr_pvalue)), mean(unlist(  list_0.5percent_fisher_jointIDA_impExpr_pvalue)))
  fisher_jointIDA_impExpr_average = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_jointIDA_impExpr_OR_average,fisher_jointIDA_impExpr_pvalue_average)
  
  
  fisher_jointIDA_OR_average = c(mean(unlist(list_0.1percent_fisher_jointIDA_OR)), mean(unlist(  list_0.2percent_fisher_jointIDA_OR )), mean(unlist(  list_0.3percent_fisher_jointIDA_OR)), mean(unlist(  list_0.4percent_fisher_jointIDA_OR)), mean(unlist(  list_0.5percent_fisher_jointIDA_OR)))
  fisher_jointIDA_pvalue_average = c(mean(unlist(  list_0.1percent_fisher_jointIDA_pvalue)), mean(unlist(  list_0.2percent_fisher_jointIDA_pvalue)), mean(unlist(  list_0.3percent_fisher_jointIDA_pvalue)), mean(unlist(  list_0.4percent_fisher_jointIDA_pvalue)), mean(unlist(  list_0.5percent_fisher_jointIDA_pvalue)))
  fisher_jointIDA_average = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_jointIDA_OR_average,fisher_jointIDA_pvalue_average)
  
  fisher_result_average = rbind(c("fisher_IDA_impExpr_average","",""),fisher_IDA_impExpr_average,c("fisher_IDA_average","",""),fisher_IDA_average,c("fisher_jointIDA_impExpr_average","",""),fisher_jointIDA_impExpr_average,c("fisher_jointIDA_average","",""),fisher_jointIDA_average)
  return(fisher_result_average)
  
}



#add by fyn 18/2/2020 caculate the median of each element in input
Fisher_Results_Median <- function(list_fisher_result)
{
  
  # list fisher_IDA_impExpr
  list_0.1percent_OR = list()
  list_0.2percent_OR = list()
  list_0.3percent_OR = list()
  list_0.4percent_OR = list()
  list_0.5percent_OR = list()
  
  list_0.1percent_pvalue = list()
  list_0.2percent_pvalue = list()
  list_0.3percent_pvalue = list()
  list_0.4percent_pvalue = list()
  list_0.5percent_pvalue = list()
  
  
  # list fisher_IDA
  list_0.1percent_fisher_IDA_OR = list()
  list_0.2percent_fisher_IDA_OR = list()
  list_0.3percent_fisher_IDA_OR = list()
  list_0.4percent_fisher_IDA_OR = list()
  list_0.5percent_fisher_IDA_OR = list()
  
  list_0.1percent_fisher_IDA_pvalue = list()
  list_0.2percent_fisher_IDA_pvalue = list()
  list_0.3percent_fisher_IDA_pvalue = list()
  list_0.4percent_fisher_IDA_pvalue = list()
  list_0.5percent_fisher_IDA_pvalue = list()
  
  
  
  # list fisher_jointIDA_impExpr
  list_0.1percent_fisher_jointIDA_impExpr_OR = list()
  list_0.2percent_fisher_jointIDA_impExpr_OR = list()
  list_0.3percent_fisher_jointIDA_impExpr_OR = list()
  list_0.4percent_fisher_jointIDA_impExpr_OR = list()
  list_0.5percent_fisher_jointIDA_impExpr_OR = list()
  
  list_0.1percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.2percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.3percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.4percent_fisher_jointIDA_impExpr_pvalue = list()
  list_0.5percent_fisher_jointIDA_impExpr_pvalue = list()
  
  
  # list fisher_jointIDA
  list_0.1percent_fisher_jointIDA_OR = list()
  list_0.2percent_fisher_jointIDA_OR = list()
  list_0.3percent_fisher_jointIDA_OR = list()
  list_0.4percent_fisher_jointIDA_OR = list()
  list_0.5percent_fisher_jointIDA_OR = list()
  
  list_0.1percent_fisher_jointIDA_pvalue = list()
  list_0.2percent_fisher_jointIDA_pvalue = list()
  list_0.3percent_fisher_jointIDA_pvalue = list()
  list_0.4percent_fisher_jointIDA_pvalue = list()
  list_0.5percent_fisher_jointIDA_pvalue = list()
  
  
  for(i in 1:length(list_fisher_result))
  {
    fisher_result = list_fisher_result[[i]]
    
    # fisher_IDA_impExpr
    fisher_IDA_impExpr = fisher_result[2:6,]
    list_0.1percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[1,2])
    list_0.2percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[2,2])
    list_0.3percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[3,2])
    list_0.4percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[4,2])
    list_0.5percent_OR[[i]] = as.numeric(fisher_IDA_impExpr[5,2])
    list_0.1percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[1,3])
    list_0.2percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[2,3])
    list_0.3percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[3,3])
    list_0.4percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[4,3])
    list_0.5percent_pvalue[[i]] = as.numeric(fisher_IDA_impExpr[5,3])
    
    
    # fisher_IDA
    fisher_IDA = fisher_result[8:12,]
    list_0.1percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[1,2])
    list_0.2percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[2,2])
    list_0.3percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[3,2])
    list_0.4percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[4,2])
    list_0.5percent_fisher_IDA_OR[[i]] = as.numeric(fisher_IDA[5,2])
    list_0.1percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[1,3])
    list_0.2percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[2,3])
    list_0.3percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[3,3])
    list_0.4percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[4,3])
    list_0.5percent_fisher_IDA_pvalue[[i]] = as.numeric(fisher_IDA[5,3])
    
    
    # fisher_jointIDA_impExpr
    fisher_jointIDA_impExpr = fisher_result[14:18,]
    list_0.1percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[1,2])
    list_0.2percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[2,2])
    list_0.3percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[3,2])
    list_0.4percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[4,2])
    list_0.5percent_fisher_jointIDA_impExpr_OR[[i]] = as.numeric(fisher_jointIDA_impExpr[5,2])
    list_0.1percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[1,3])
    list_0.2percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[2,3])
    list_0.3percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[3,3])
    list_0.4percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[4,3])
    list_0.5percent_fisher_jointIDA_impExpr_pvalue[[i]] = as.numeric(fisher_jointIDA_impExpr[5,3])
    
    
    # fisher_jointIDA
    fisher_jointIDA = fisher_result[20:24,]
    list_0.1percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[1,2])
    list_0.2percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[2,2])
    list_0.3percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[3,2])
    list_0.4percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[4,2])
    list_0.5percent_fisher_jointIDA_OR[[i]] = as.numeric(fisher_jointIDA[5,2])
    list_0.1percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[1,3])
    list_0.2percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[2,3])
    list_0.3percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[3,3])
    list_0.4percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[4,3])
    list_0.5percent_fisher_jointIDA_pvalue[[i]] = as.numeric(fisher_jointIDA[5,3])
  }
  
  
  
  fisher_IDA_impExpr_OR_median = c(median(unlist(list_0.1percent_OR)),median(unlist(list_0.2percent_OR)),median(unlist(  list_0.3percent_OR)),median(unlist(  list_0.4percent_OR)),median(unlist(  list_0.5percent_OR )))
  fisher_IDA_impExpr_pvalue_median = c(median(unlist(list_0.1percent_pvalue)),median(unlist(list_0.2percent_pvalue)),median(unlist(  list_0.3percent_pvalue)),median(unlist(  list_0.4percent_pvalue)),median(unlist(  list_0.5percent_pvalue)))
  
  fisher_IDA_impExpr_median = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_IDA_impExpr_OR_median,fisher_IDA_impExpr_pvalue_median)
  
  
  fisher_IDA_OR_median = c(median(unlist(list_0.1percent_fisher_IDA_OR)), median(unlist(  list_0.2percent_fisher_IDA_OR )), median(unlist(  list_0.3percent_fisher_IDA_OR)), median(unlist(  list_0.4percent_fisher_IDA_OR)), median(unlist(  list_0.5percent_fisher_IDA_OR)))
  fisher_IDA_pvalue_median = c(median(unlist(  list_0.1percent_fisher_IDA_pvalue)), median(unlist(  list_0.2percent_fisher_IDA_pvalue)), median(unlist(  list_0.3percent_fisher_IDA_pvalue)), median(unlist(  list_0.4percent_fisher_IDA_pvalue)), median(unlist(  list_0.5percent_fisher_IDA_pvalue)))
  fisher_IDA_median = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_IDA_OR_median,fisher_IDA_pvalue_median)
  
  fisher_jointIDA_impExpr_OR_median = c(median(unlist(list_0.1percent_fisher_jointIDA_impExpr_OR)), median(unlist(  list_0.2percent_fisher_jointIDA_impExpr_OR )), median(unlist(  list_0.3percent_fisher_jointIDA_impExpr_OR)), median(unlist(  list_0.4percent_fisher_jointIDA_impExpr_OR)), median(unlist(  list_0.5percent_fisher_jointIDA_impExpr_OR)))
  fisher_jointIDA_impExpr_pvalue_median = c(median(unlist(  list_0.1percent_fisher_jointIDA_impExpr_pvalue)), median(unlist(  list_0.2percent_fisher_jointIDA_impExpr_pvalue)), median(unlist(  list_0.3percent_fisher_jointIDA_impExpr_pvalue)), median(unlist(  list_0.4percent_fisher_jointIDA_impExpr_pvalue)), median(unlist(  list_0.5percent_fisher_jointIDA_impExpr_pvalue)))
  fisher_jointIDA_impExpr_median = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_jointIDA_impExpr_OR_median,fisher_jointIDA_impExpr_pvalue_median)
  
  
  fisher_jointIDA_OR_median = c(median(unlist(list_0.1percent_fisher_jointIDA_OR)), median(unlist(  list_0.2percent_fisher_jointIDA_OR )), median(unlist(  list_0.3percent_fisher_jointIDA_OR)), median(unlist(  list_0.4percent_fisher_jointIDA_OR)), median(unlist(  list_0.5percent_fisher_jointIDA_OR)))
  fisher_jointIDA_pvalue_median = c(median(unlist(  list_0.1percent_fisher_jointIDA_pvalue)), median(unlist(  list_0.2percent_fisher_jointIDA_pvalue)), median(unlist(  list_0.3percent_fisher_jointIDA_pvalue)), median(unlist(  list_0.4percent_fisher_jointIDA_pvalue)), median(unlist(  list_0.5percent_fisher_jointIDA_pvalue)))
  fisher_jointIDA_median = cbind(c(0.1,0.2,0.3,0.4,0.5),fisher_jointIDA_OR_median,fisher_jointIDA_pvalue_median)
  
  fisher_result_median = rbind(c("fisher_IDA_impExpr_median","",""),fisher_IDA_impExpr_median,c("fisher_IDA_median","",""),fisher_IDA_median,c("fisher_jointIDA_impExpr_median","",""),fisher_jointIDA_impExpr_median,c("fisher_jointIDA_median","",""),fisher_jointIDA_median)
  return(fisher_result_median)
  
}

