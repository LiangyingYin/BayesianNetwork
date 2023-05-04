
library(dplyr)
library(data.table)

get_outcome_samples_for_binary_trait <- function(outcome,target_outcome_name){
  cases = outcome[outcome[,target_outcome_name]==1,]
  ctrls = outcome[outcome[,target_outcome_name]==0,]
  #selected_ctrls_index = sample(1:nrow(ctrls),size = nrow(cases))
  #selected_ctrls = ctrls[selected_ctrls_index,]
  #results = rbind(cases,selected_ctrls)
  results = rbind(cases,ctrls)
  return(results)
}

#*****************************************************************
# Revised by: yinly
# Date: Feb 27, 2020
# Description: get samples for binary traits
#*****************************************************************
get_outcome_samples_for_Diabete <- function(outcome,target_outcome_name){
  cases = outcome[outcome[,target_outcome_name]==1,]
  ctrls = outcome[outcome[,target_outcome_name]==0,]
  #selected_ctrls_index = sample(1:nrow(ctrls),size = nrow(cases))
  #selected_ctrls = ctrls[selected_ctrls_index,]
  #results = rbind(cases,selected_ctrls)
  results = rbind(cases,ctrls)
  return(results)
}
#*************************************************************************
# Perform univariate test
# Parameters:
# target_outcome_name: index of studied trait in UKBB phenotype file
# trait: studied trait(corresponding with target_outcome_name)
# tissue: causal tissue for the studied trait
#*************************************************************************
get_univariate_test_results <-function(trait,tissue,target_outcome_name,outcome_index){
  #target_outcome_name = "CAD"
  expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Whole_Blood_CIS_Trans_expression.txt"
  #expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Lung_predicted_expression.txt"
  pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_31March2022.txt"
  #___  END OF INPUT__________________________________________
  
  expr = fread(expression_file_dir )
  
  #************************************
  # screen away genes with very low variance (this will lead to NA entries after scaling in huge
  #***********************************
  #var_thres = 0.01
  #var_list = apply(expr[,3:ncol(expr)], 2, var)
  #gene_expr = data.frame(expr[,3:ncol(expr)])
  #gene_expr_filtered =  gene_expr[,var_list>var_thres]
  #expr <- cbind(expr[,1:2],gene_expr_filtered )
  
  #********************************
  # Extract outcome trait - Just take BMI as an example; code can be searched from the HTML file
  #********************************
  pheno = fread(pheno_dir)
  outcome = cbind(pheno[,1], pheno[,..target_outcome_name])
  outcome  = data.frame(   na.omit(outcome)   )
  
  print(dim(outcome))
  
  #***********************************************************************
  # Revised by: yinly
  # Date: Feb 27, 2020
  # Description: select cases and controls for studied traits for diabetes
  #***********************************************************************
  if(trait =="Diabetes"){
    #outcome_index = trait
    outcome = get_outcome_samples_for_Diabete(outcome,outcome_index)
  }
  #outcome = get_outcome_samples_for_binary_trait(outcome,outcome_index)
  
  #**************************
  #                 merging expression with the outcome
  #******************************
  expr = data.frame(expr)
  # Here, after apply the inner_join function, the column name of outcome changed to "X30690.0.0"(original:30690-0.0)
  df = inner_join(outcome, expr, by="FID")   #merge(outcome, expr, by="FID")
  print(  (dim(df))  )
  
  
  #*****************************************
  # first regress y and x on the top 10 PCs respectively to correct for population stratification
  #******************************************
  #extract principal components from UKBB phenotype file
  # PC_start_ind = which( colnames(pheno) == "X22009.0.1")
  # Revised by yinly on Jan 6, 2020
  # Based on the newly extracted phenotypes file, the column name for the 1st PC has changed to 22009.0.1
  
  pheno = data.frame(pheno)
  PC_start_ind = which( colnames(pheno) == "X22009.0.1")
  PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
  
  # PCs_prefix = "22009.0."
  # PCs_suffix = seq(1,10)
  # PC_names = paste0(PCs_prefix,PCs_suffix)
  # PCs = cbind(pheno[,"FID"], pheno[,PC_names])
  colnames(PCs)[1] <- "FID"
  colnames(PCs)[2:11] <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")
  
  #***********************************************************************
  # Revised by: yinly, Date: March 6, 2020
  # Description: preselect a set of genes that are related to the outcome
  #***********************************************************************
  df = inner_join(df,PCs,by ="FID")
  #df = univariate_gene_selection(df,p_thres,outcome_index)
  #print(df)
  #df = inner_join(df, PCs, by ="FID")
  
  
  
  #*****************
  # count the no. of genes available
  #******************
  grep_res = grep("ENSG", colnames(df) )
  start.ind = min(grep_res)
  end.ind = max(grep_res)
  no_genes = end.ind - start.ind + 1
  ##start.ind = index for which we start to have gene expressions
  binary = TRUE
  univariate.mat = matrix(nrow = no_genes, ncol= 3 )
  attach(df)
  # Residualize y on PCs
  lm.obj_y = lm(df[,outcome_index] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
  df[,outcome_index] = lm.obj_y$resid
  for (i in 1:no_genes) {
    # Residualize on each gene
    lm.obj_x = lm(df[,start.ind+ i-1] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
    df[,start.ind+ i-1] = lm.obj_x$resid
    # Commented by:yinly
    # Perform univariate test by linear regression to select significant genes
    #lm.obj = lm(df[,outcome_index] ~ df[,start.ind+ i-1] + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
    lm.obj = lm(df[,outcome_index] ~ df[,start.ind+ i-1])
    #resid.mat[,i] = lm.obj$resid
    univariate.mat[i,1] = colnames(df)[start.ind + i - 1]
    univariate.mat[i,2] = summary(lm.obj)$coeff[2,1]
    univariate.mat[i,3] = summary(lm.obj)$coeff[2,4]
  }
  
  colnames(univariate.mat) = c("Genes","Estimates","Pvalues")
  #tissue = "Adipose_Subcutaneous"
  Univariate_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Univariate_Results_All_Rerun/"
  Univariate_filepath = paste0(Univariate_prefix,tissue,"_",trait,"_Univariate_Linear.csv")
  fwrite(as.data.frame(univariate.mat),Univariate_filepath)
}

# Load Dictionary for studied traits and their causal tissues
Dict_filepath = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Whole_Blood_associated_traits.csv"
Dict = as.matrix(fread(Dict_filepath))
#for(i in 1:10){
  #i = 21
  trait = "COVID_Confirm"
  tissue = "Whole_Blood"
  target_outcome_name = "COVID_Confirm"
  outcome_index = "COVID_Confirm"
  #trait = Dict[i,1]
  #tissue = Dict[i,2]
  #target_outcome_name = Dict[i,3]
  #outcome_index = Dict[i,4]
  #cat("Univariate test for the",i,"trait",trait,"with",tissue,"as causal tissue","\n")
  cat("Univariate test for the trait",trait,"with",tissue,"as causal tissue","\n")
  get_univariate_test_results(trait,tissue,target_outcome_name,outcome_index)
  cat("Univariate test for",trait,"with",tissue,"as causal tissue is completed!","\n")
#}

