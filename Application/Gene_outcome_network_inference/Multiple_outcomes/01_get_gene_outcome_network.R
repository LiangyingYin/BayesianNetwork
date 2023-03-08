
#source("/exeh/exe3/yinly/BayesianNetwork/03_UKBB/src/PCselect_rawData_stable.R")
source("/exeh_3/yinly//BayesianNetwork/03_UKBB_Network/src/PCSelect_Parallel.R")
library(dplyr)
library(data.table)
library(RNOmni)
library(parallel)
#library(biglm)

#*****************************************************************
# Revised by: yinly
# Date: Feb 27, 2020
# Description: get samples for binary traits
#*****************************************************************
get_outcome_samples_for_binary_trait <- function(outcome,target_outcome_name){
  cases = outcome[outcome[,target_outcome_name]==1,]
  ctrls = outcome[outcome[,target_outcome_name]==0,]
  #selected_ctrls_index = sample(1:nrow(ctrls),size = nrow(cases))
  #selected_ctrls = ctrls[selected_ctrls_index,]
  #results = rbind(cases,selected_ctrls)
  results = rbind(cases,ctrls)
  return(results)
}


univariate_gene_selection <- function(df,p_thres,outcome_index){
  #*****************
  # count the no. of genes available
  #******************
  grep_res = grep("ENSG", colnames(df) )
  start.ind = min(grep_res)
  end.ind = max(grep_res)
  no_genes = end.ind - start.ind + 1
  ##start.ind = index for which we start to have gene expressions
  genes_assoc = matrix(nrow = no_genes,ncol=2) 
  attach(df)
  p1=proc.time()
  for (i in 1:no_genes){
    glm.obj = glm(df[,outcome_index] ~ df[,start.ind+ i-1]+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial)
    #glm.obj <- bigglm(df[,outcome_index] ~ df[,start.ind+ i-1]+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial(), chunksize=100, maxit=10)
    genes_assoc[i,1] = start.ind+ i-1
    genes_assoc[i,2] = summary(glm.obj)$coeff[2,4]
  }
  proc.time()-p1
  #summary(genes_assoc)
  sig_genes_index = genes_assoc[(genes_assoc[,2]<p_thres),1]
  detach(df)
  FID_index = grep("FID", colnames(df) )
  PC_index = grep("PC",colnames(df))
  Trait_index = grep(outcome_index,colnames(df))
  selected_feats_index =  c(FID_index,sig_genes_index,PC_index,Trait_index)
  df = df[,selected_feats_index]
  return(df)
}

setwd("/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/02_Multiple_Traits/Order3_Results_CIS_Trans/")
get_Causal_Network_Based_On_UKBB <- function(target_outcome_name,outcome_index,outcome_name,p_thres){
  
  # Input the coded ID for interested variable
  # target_outcome_name = "X21001.0.0"  #Just take BMI as an example; code can be searched from the HTML file
  # target_outcome_name = "X20127.0.0" # Basic metabolic rate
  # Revised by yinly on Jan 6,2020
  # For details about the column names of the phenotypes, please refer to:
  # UKBB_Phenotypes_Dict_for_Extracted_Continuous_and_Binary_Traits.xlsx
  # target_outcome_name = "30690-0.0" # This is the column name for cholestral
  # pheno_dir = "/exeh/exe4/sohc/UKBB2/data/phenotype/array_phenos.txt"
  # Revised by yinly on Jan 6, 2020
  # Use new extracted phenotypes
  var_thres = 0.01
  expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/UKBB_Whole_Blood_CIS_Trans_expression.txt"
  predicted_pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/02_Multiple_Traits/Genetically_Predicted_BMI.csv"
  #pheno_dir = "/exeh_3/yinly/UKBB/01_Data/UKBB_Traits_with_ID_July29.csv"
  #pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_with_FID_updated.csv"
  pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Cov19.txt"
  #___  END OF INPUT__________________________________________
  
  expr = fread(expression_file_dir )
  #************************************
  # screen away genes with very low variance (this will lead to NA entries after scaling in huge
  #***********************************
  var_list = apply(expr[,3:ncol(expr)], 2, var)
  gene_expr = data.frame(expr[,3:ncol(expr)])
  gene_expr_filtered =  gene_expr[,var_list>var_thres]
  expr <- cbind(expr[,1:2],gene_expr_filtered )

  #************************************
  # Load genetically predicted pheno and combine with gene expression data
  #************************************
  predicted_pheno = fread(predicted_pheno_dir)
  colnames(predicted_pheno) = c("FID","ENSG_BMI")
  expr_new = inner_join(expr, predicted_pheno,by="FID")
  
  #********************************
  # Extract outcome trait - Just take BMI as an example; code can be searched from the HTML file
  #********************************
  pheno = fread(pheno_dir)
  outcome = cbind(pheno[,1], pheno[,..target_outcome_name])
  outcome  = data.frame(   na.omit(outcome)   )
  #*****************************************************************
  # Revised by: yinly
  # Date: Feb 27, 2020
  # Description: select controls based on the sample size of cases
  #*****************************************************************
  outcome = get_outcome_samples_for_binary_trait(outcome,outcome_index)
  #pheno_brief = fread("/exeh/exe4/sohc/UKBB2/data/phenotype/UKBBphenos.txt")   ##formatted file by Timothy Mak
  
  # male_data = fread("/exeh_3/yinly/Cardiovascular/Data/Male_Finland_RFImputed_Phenotypes.csv")
  # fem_data = fread("/exeh_3/yinly/Cardiovascular/Data/Female_Finland_RFImputed_Phenotypes.csv")
  # pheno = rbind(male_data, fem_data)
  
  
  
  
  #**************************
  #                 merging expression with the outcome
  #******************************
  expr = data.frame(expr_new)
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
  df = univariate_gene_selection(df,p_thres,outcome_index)
  
  #df = inner_join(df, PCs, by ="FID")
  
  
  
  #*****************
  # count the no. of genes available
  #******************
  grep_res = grep("ENSG", colnames(df) )
  start.ind = min(grep_res)
  end.ind = max(grep_res)
  no_genes = end.ind - start.ind + 1
  ##start.ind = index for which we start to have gene expressions
  
  #********************************
  # Residualized X
  #*********************************
  resid.mat = matrix(nrow = nrow(df), ncol= no_genes )
  attach(df)
  for (i in 1:no_genes) {
    lm.obj = lm(df[,start.ind+ i-1] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
    resid.mat[,i] = lm.obj$resid
  }
  
  colnames(resid.mat) <- colnames(df)[start.ind:end.ind]
  
  #********************************
  # Residualized Y
  #*********************************
  # Revised by yinly on Jan 7,2020,the inner_join function changed the column name of the target_outcome_name
  # # Reset the value of
  # outcome_index = "X30690.0.0"
  lm.obj_y = lm(df[,outcome_index] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
  resid_y = lm.obj_y$resid
  detach(df)
  
  dfnew = data.frame(outcome= resid_y,  resid.mat)
  #save(dfnew, file= "df_new_original_wholeBlood_UKBB.Rdata")
  #*************************
  # Rank-based Inverse normal transformaton
  #************************
  dfnew_rankNorm = apply(dfnew, 2, rankNorm)
  #save(dfnew_rankNorm, file= "df_new_rankNorm_wholeBlood_UKBB.Rdata")
  
  
  
  #**********************************
  # Pre-computing correlation matrix
  # We shall just compute the cor mat for gene expressions 1st, as it can be reused for different outcomes
  #*******************************
  library(coop)
  #load("df_new_rankNorm_wholeBlood_UKBB.Rdata")
  
  #**************************************
  #  IMPORTANT: the following expression corr matrix can be re-used for other outcomes (given we are working on the SAME tissue)
  #**************************************
  #******************************************
  # Revised by: yinly, Date: March 6, 2020
  #******************************************
  expr_corMat = pcor(dfnew_rankNorm[,-1])
  # save(expr_corMat, file="correl_mat_expr_rankNorm_UKBB_wholeBlood.Rdata")
  #************************************************************************
  # Since we have calculated the corrlation matrix for whole blood tissue,
  # we could just directly reload the calculated matrix
  #***********************************************************************
  #load("/exeh/exe4/sohc/BN_to_UKBB/correl_mat_expr_rankNorm_UKBB_wholeBlood.Rdata")
  
  #now add back the correlation between the outcome and each feature
  #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
  cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
  precompute_corMat = cbind(cor_with_outcome, expr_corMat)
  
  precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat )
  
  
  #********************************************************
  # Apply PC-simple algorithm around the response variable on NFBC
  #**************************************************************
  library(pcalg)
  #library(kpcalg)
  
  ## Load  data
  n <- nrow (dfnew_rankNorm)
  V <- colnames(dfnew_rankNorm) # labels aka node names
  
  t1 = proc.time()
  ## estimate local causal network structure around the response variable
  
  # Revise by yinly:
  # Date: Feb 12,2020
  # Description: set more stringent alpha for binary traits
  pcSimple.fit <- PCSelect_Parallel(y=dfnew_rankNorm[,1],
                                    dm=dfnew_rankNorm[,-1] ,
                                    method = c("parallel"),
                                    mem.efficient = FALSE,
                                    num_workers = 10,
                                    #alpha=0.001,
                                    alpha=0.001,
                                    corMat = precompute_corMat,
                                    max_ord=3,
                                    corMethod = "standard",
                                    verbose = TRUE, directed = TRUE)
  
  proc.time()-t1
  filepath = paste0("UKBB_BMI_",outcome_name,"_Whole_Blood.Rdata")
  save(pcSimple.fit, file = filepath)
  #save(pcSimple.fit, file="UKBB_Creatine.Rdata")
}

library(data.table)

target_outcome_name = "2443-0.0" #2443-0.0 for Diabetes # 30750-0.0 for HbA1c
outcome_index = "X2443.0.0"
outcome_name = "Diabetes"
cat("This is the trait:",outcome_name,"\n")
get_Causal_Network_Based_On_UKBB(target_outcome_name,outcome_index,outcome_name,0.2)
cat("The trait:",outcome_name,"is completed!\n")
#Traits = as.matrix(fread("/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/Binary_Traits_Dict.csv",header = TRUE))
# Define studied trait
# CAD
#i = 6
#for(i in 1:nrow(Traits)){
#  i = 6
#  cat("This is the",i,"trait:",Traits[i,3],"\n")
#  get_Causal_Network_Based_On_UKBB(Traits[i,1],Traits[i,2],Traits[i,3],0.2)
#  cat("The",i,"trait:",Traits[i,3],"is completed!","\n")
#}
