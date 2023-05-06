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
# to do list: Code optimization
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

# #fisher_result_single is a matrix 5*3
# Fisher_Results_Average_subfunction <- function(fisher_result_single)
# {
#   
# }


#######################temp test for using####################################
Get_dtIDA <- function(tureGraph,cov_imputedExpr,result_list,pcSimple.fit,genes_num,pc.fit)
{
  #try to convert graphNEL to igraph in order to obtain the edgelist.
  
  nodes_tureGraph = nodes(tureGraph)
  phenotype_loc = length(nodes_tureGraph)  # phenotype is the last node in the Graph
  nodes(tureGraph) = as.character(seq(1,phenotype_loc))
  nodes_tureGraph = nodes(tureGraph)
  
  dt_IDA = data.frame()
  max_IDAlength1 = 0
  
  
  
  
  for(gene in 1:genes_num)
  {
    node1 = nodes_tureGraph[gene] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    
    covTrue <- trueCov(tureGraph)
    true_IDA <- ida(node1,node2, covTrue, tureGraph, method = "local", type = "pdag")
    
    l.IDA.estimate <- ida(node1,node2, cov_imputedExpr , pc.fit, method = "local", type = "pdag")									   
    
    
    print(l.IDA.estimate) # print the estimate IDA value
    
    num_local_IDA = length(l.IDA.estimate)
    if(num_local_IDA > max_IDAlength1)
    {
      max_IDAlength1 = num_local_IDA
    }
  }
  
  max_IDAlength = max_IDAlength1
  #tureGraph as standard
  
  
  nodes_true = nodes(tureGraph)
  x.pos = nodes_true[!nodes_true %in% node2]
  y.pos = node2
  
  x.pos = as.numeric(x.pos) # change the data type from character to numeric
  y.pos = as.numeric(y.pos)
  
  #jointIDA compute
  true_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,covTrue,graphEst=tureGraph,technique="RRC")
  
  #estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov(dat2),graphEst=estimateGraph,technique="RRC")
  estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov_imputedExpr,graphEst= pc.fit,technique="RRC")
  
  
  for(gene in 1:genes_num)
  {
    node1 = nodes_tureGraph[gene] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    
    #IDA compute
    covTrue <- trueCov(tureGraph)
    true_IDA <- ida(node1,node2, covTrue, tureGraph, method = "local", type = "pdag")
    # true_IDA = causalEffect(tureGraph,node2,node1)
    
    
    l.IDA.estimate <- ida(node1,node2, cov_imputedExpr , pc.fit, method = "local", type = "pdag")	
    
    jointIDA_position = which(x.pos==node1)
    true_jointIDA = true_jointIDA_list[jointIDA_position]
    
    
    dt_IDA[gene,1] = node1 # first node in the edge
    dt_IDA[gene,2] = node2 # second node in the edge
    
    dt_IDA[gene,3] = true_IDA
    num_local_IDA = length(l.IDA.estimate)
    
    
    for(numlIDA in 1:num_local_IDA)
    {
      dt_IDA[gene,(3+numlIDA)] = l.IDA.estimate[numlIDA]
    }
    
    #filled missing value
    tmp1=  4+max_IDAlength
    tmp2 = 3+numlIDA
    minus_1_2 = tmp1-tmp2
    if(minus_1_2>1)
    {
      for(tmp in 1:minus_1_2)
      {
        # dt_IDA[gene_num,(tmp2+tmp-1)] = 'NA'
        dt_IDA[gene,(tmp2+tmp)] = 'NA'
      }
    }
    
    
    dt_IDA[gene,(4+max_IDAlength)] = true_jointIDA
    
    num_estimate_jointIDA = ncol(estimate_jointIDA_list)
    for(num_jointIDA in 1:num_estimate_jointIDA)
    {
      dt_IDA[gene,(4+max_IDAlength+num_jointIDA)] = estimate_jointIDA_list[jointIDA_position,num_jointIDA]
    }
  }
  
  colnames(dt_IDA)[1] = "gene_loc"
  colnames(dt_IDA)[2] = "phenotype_loc"
  colnames(dt_IDA)[3] = "true_IDA"
  colnames(dt_IDA)[4] = "Estimate_IDA"
  colnames(dt_IDA)[(4+max_IDAlength)] = "true_jointIDA"
  
  
  print("This is dt_IDA")
  print(dt_IDA)
}

Get_dtIDA2 <- function(tureGraph,cov_rowExpr,cov_imputedExpr,result_list,pcSimple.fit,genes_num,pc.fit)
{
  #try to convert graphNEL to igraph in order to obtain the edgelist.
  nodes_tureGraph = nodes(tureGraph)
  phenotype_loc = length(nodes_tureGraph)  # phenotype is the last node in the Graph
  nodes(tureGraph) = as.character(seq(1,phenotype_loc))
  nodes_tureGraph = nodes(tureGraph)
  
  dt_IDA = data.frame()
  max_IDAlength1 = 0
  
  for(gene in 1:genes_num)
  {
    node1 = nodes_tureGraph[gene] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    
    covTrue <- trueCov(tureGraph)
    true_IDA <- ida(node1,node2, covTrue, tureGraph, method = "local", type = "pdag")
    
    l.IDA.estimate <- ida(node1,node2, cov_imputedExpr , pc.fit, method = "local", type = "pdag")							   
    
    
    print(l.IDA.estimate) # print the estimate IDA value
    
    num_local_IDA = length(l.IDA.estimate)
    if(num_local_IDA > max_IDAlength1)
    {
      max_IDAlength1 = num_local_IDA
    }
  }
  
  max_IDAlength = max_IDAlength1
  #tureGraph as standard
  
  
  nodes_true = nodes(tureGraph)
  x.pos = nodes_true[!nodes_true %in% node2]
  y.pos = node2
  
  x.pos = as.numeric(x.pos) # change the data type from character to numeric
  y.pos = as.numeric(y.pos)
  
  #jointIDA compute
  true_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,covTrue,covTrue,graphEst=tureGraph,technique="RRC")
  
  #estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov(dat2),graphEst=estimateGraph,technique="RRC")
  estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov_rowExpr,cov_imputedExpr,graphEst = pc.fit,technique="RRC")
  
  
  for(gene in 1:genes_num)
  {
    node1 = nodes_tureGraph[gene] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    
    #IDA compute
    covTrue <- trueCov(tureGraph)
    true_IDA <- ida(node1,node2, covTrue, tureGraph, method = "local", type = "pdag")
    # true_IDA = causalEffect(tureGraph,node2,node1)
    
    
    l.IDA.estimate <- ida(node1,node2, cov_imputedExpr , pc.fit, method = "local", type = "pdag")	
    
    jointIDA_position = which(x.pos==node1)
    true_jointIDA = true_jointIDA_list[jointIDA_position]
    
    
    dt_IDA[gene,1] = node1 # first node in the edge
    dt_IDA[gene,2] = node2 # second node in the edge
    
    dt_IDA[gene,3] = true_IDA
    num_local_IDA = length(l.IDA.estimate)
    
    
    for(numlIDA in 1:num_local_IDA)
    {
      dt_IDA[gene,(3+numlIDA)] = l.IDA.estimate[numlIDA]
    }
    
    #filled missing value
    tmp1=  4+max_IDAlength
    tmp2 = 3+numlIDA
    minus_1_2 = tmp1-tmp2
    if(minus_1_2>1)
    {
      for(tmp in 1:minus_1_2)
      {
        # dt_IDA[gene_num,(tmp2+tmp-1)] = 'NA'
        dt_IDA[gene,(tmp2+tmp)] = 'NA'
      }
    }
    
    
    dt_IDA[gene,(4+max_IDAlength)] = true_jointIDA
    
    num_estimate_jointIDA = ncol(estimate_jointIDA_list)
    for(num_jointIDA in 1:num_estimate_jointIDA)
    {
      dt_IDA[gene,(4+max_IDAlength+num_jointIDA)] = estimate_jointIDA_list[jointIDA_position,num_jointIDA]
    }
  }
  
  colnames(dt_IDA)[1] = "gene_loc"
  colnames(dt_IDA)[2] = "phenotype_loc"
  colnames(dt_IDA)[3] = "true_IDA"
  colnames(dt_IDA)[4] = "Estimate_IDA"
  colnames(dt_IDA)[(4+max_IDAlength)] = "true_jointIDA"
  
  
  print("This is dt_IDA")
  print(dt_IDA)
}


CVP_rho_matrix2 <- function(X = NULL,
                            rho_matrix ,
                            lam = 10^seq(-2, 1, 0.1),
                            diagonal = FALSE,
                            tol = 1e-04,
                            maxit = 10000,
                            adjmaxit = NULL,
                            K = 5,
                            crit.cv = c("loglik",
                                        "AIC", "BIC", "eBIC"),
                            ebic_gamma = 0.5,
                            start = c("warm", "cold"),
                            cores = 5, trace = c("progress",
                                                 "print", "none"), ...) {
  
  # match values
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  
  # make cluster and register cluster
  num_cores = detectCores()
  if (cores > num_cores) {
    cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
  }
  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }
  
  cluster = makeCluster(cores)
  registerDoParallel(cluster)
  
  # use cluster for each fold in CV
  n = nrow(X)
  ind = sample(n)
  k = NULL
  CV = foreach(k = 1:K, .packages = c("CVglasso","glassoFast"), .combine = "cbind",
               .inorder = FALSE) %dopar% {
                 
                 # set progress bar
                 if (trace == "progress") {
                   progress = txtProgressBar(max = length(lam), style = 3)
                 }
                 
                 # training set
                 leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * n/K)]
                 X.train = X[-leave.out, , drop = FALSE]
                 X_bar = apply(X.train, 2, mean)
                 X.train = scale(X.train, center = X_bar, scale = FALSE)
                 
                 # validation set
                 X.valid = X[leave.out, , drop = FALSE]
                 X.valid = scale(X.valid, center = X_bar, scale = FALSE)
                 
                 # sample covariances
                 S.train = crossprod(X.train)/(dim(X.train)[1])
                 S.valid = crossprod(X.valid)/(dim(X.valid)[1])
                 
                 # initial estimates
                 init = S.train
                 initOmega = diag(ncol(S.train))
                 CV_error = array(0, length(lam))
                 
                 # initial sigma
                 if (!diagonal) {
                   
                   # provide estimate that is pd and dual feasible
                   Sminus = S.train
                   diag(Sminus) = 0
                   alpha = min(c(lam[1]/max(abs(Sminus)), 1))
                   init = (1 - alpha) * S.train
                   diag(init) = diag(S.train)
                   
                 }
                 
                 # loop over all tuning parameters
                 for (i in 1:length(lam)) {
                   
                   # set temporary tuning parameter
                   lam_ = lam[i]
                   
                   # update diagonal elements of init, if necessary
                   if (diagonal) {
                     diag(init) = diag(S.train) + lam_
                   }
                   
                   # compute the penalized likelihood precision matrix estimator
                   #_________________________________________________________
                   # HC: changed rho = lam_ to 'rho = rho_matrix*lam_'
                   #________________________________________________________
                   # GLASSO = glasso(s = S.train, rho = rho_matrix*lam_ , thr = tol, maxit = maxit,
                   #     penalize.diagonal = diagonal, start = "warm", w.init = init,
                   #     wi.init = initOmega, trace = FALSE)
                   GLASSO = glassoFast(S = S.train, rho = rho_matrix*lam_ , thr = tol,
                                       start = "warm", w.init = init, wi.init = initOmega,
                                       trace = FALSE)
                   if (start == "warm") {
                     
                     # option to save initial values for warm starts
                     init = GLASSO$w
                     initOmega = GLASSO$wi
                     maxit = adjmaxit
                     
                   }
                   
                   # compute the observed negative validation loglikelihood (close
                   # enoug)
                   
                   # Cf iDINGO R package https://rdrr.io/cran/iDINGO/src/R/extendedBIC.R
                   
                   
                   CV_error[i] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) -
                                                  determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
                   #if we use GLASSO$w instead of GLASSO$wi, it seems work, we can get the min lamada, which in crease cor(trueIDA&estimateIDA)
                   # CV_error[i] = (nrow(X)/2) * (sum(GLASSO$w * S.valid) -
                   #                                determinant(GLASSO$w, logarithm = TRUE)$modulus[1])
                   
                   # CV_error[i] = CV_error[i] + sum(GLASSO$w != 0)
                   
                   #_____________________________________________
                   ##amended in May 2020 
                   # https://arxiv.org/pdf/1011.6640.pdf;  http://www3.stat.sinica.edu.tw/sstest/oldpdf/A22n310.pdf
                   #no. of edges used isntead of sum(GLASSO$wi) 
                   #_________________________________________________
                   
                   # update for crit.cv, if necessary
                   if (crit.cv == "AIC") {
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0)
                     # CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0)
                   }
                   
                   # this is the traditiona BIC/2
                   # http://www3.stat.sinica.edu.tw/sstest/oldpdf/A22n310.pdf
                   if (crit.cv == "BIC") {
                     # CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0) *
                     #   log(nrow(X))/2
                     
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0) *
                       log(nrow(X))/2
                   }
                   
                   #_________________________________
                   # added May 2020
                   #__________________________________
                   # this is the extended BIC/2  #https://arxiv.org/pdf/1011.6640.pdf ;
                   if (crit.cv == "eBIC") {
                     # CV_error[i] = CV_error[i] + sum(GLASSO$wi != 0) *
                     #   log(nrow(X))/2
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0) *
                       log(nrow(X))/2 + 2*sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0)* ebic_gamma * log(ncol(X))
                   }
                   
                   
                   # update progress bar
                   if (trace == "progress") {
                     setTxtProgressBar(progress, i)
                     
                     # if not quiet, then print progress lambda
                   } else if (trace == "print") {
                     cat("\nFinished lam = ", paste(lam_, sep = ""))
                   }
                 }
                 
                 # return CV errors
                 return(CV_error)
                 
               }
  
  # determine optimal tuning parameters
  AVG = apply(CV, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  
  # stop cluster
  stopCluster(cluster)
  
  # return best lam and alpha values
  return(list(lam = best_lam, min.error = error, avg.error = AVG,
              cv.error = CV))
  
}