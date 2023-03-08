
#source("/exeh/exe3/yinly/BayesianNetwork/03_UKBB/src/PCselect_rawData_stable.R")
source("/exeh_3/yinly//BayesianNetwork/03_UKBB_Network/src/PCSelect_Parallel.R")
library(dplyr)
library(data.table)
library(RNOmni)
library(parallel)
#library(biglm)

#*****************************************************************
# Revised by: yinly
# Date: July 9, 2020
# Description: get samples for cov19
# Revised by yinly: August 26,2020
# Description: utilzing subsampling strategy on the dataset
#*****************************************************************

outcome_subsampling_for_binary_trait <-function(outcome,target_outcome_name,sampling_ratio){
  ctrls = outcome[outcome[,target_outcome_name]==0,]
  size_sampled_ctrls = round(nrow(ctrls)*sampling_ratio)
  sampled_ctrls_index = sample(1:nrow(ctrls),size = size_sampled_ctrls)
  sampled_ctrls = ctrls[sampled_ctrls_index,]
  cases = outcome[outcome[,target_outcome_name]==1,]
  size_sampled_cases = round(nrow(cases)*sampling_ratio)
  sampled_cases_index = sample(1:nrow(cases),size = size_sampled_cases)
  sampled_cases = cases[sampled_cases_index,]
  results = rbind(sampled_cases,sampled_ctrls)
  return(results)
}
outcome_subsampling_for_continuous_trait <-function(outcome,sampling_ratio){
  size_sampled = round(nrow(outcome)*sampling_ratio)
  sampled_index = sample(1:nrow(outcome),size = size_sampled)
  sampled_outcome = outcome[sampled_index,]
  return(sampled_outcome)
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
    #glm.obj = glm(df[,outcome_index] ~ df[,start.ind+ i-1]+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial)
    #glm.obj <- bigglm(df[,outcome_index] ~ df[,start.ind+ i-1]+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial(), chunksize=100, maxit=10)
    lm.obj = lm(df[,outcome_index] ~ df[,start.ind+ i-1]+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
    genes_assoc[i,1] = start.ind+ i-1
    #genes_assoc[i,2] = summary(glm.obj)$coeff[2,4]
    genes_assoc[i,2] = summary(lm.obj)$coeff[2,4]
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

Get_Bootstrapped_Gene_Outcome_Network <- function(expr,outcome,target_outcome_name,sampling_ratio,iter,PCs,sampling_cutoff,p_thres){
  # Define assoc_mat and zMin_mat to store the results for each bootstrap
  assoc_mat = NULL
  zMin_mat = NULL
  qlambda_sum = 0 # This is used for the calculation of E(V)
  PCSelect_Result_final_list = list()
  for(i in 1:iter){
    #outcome = outcome_subsampling_for_cov19(outcome,target_outcome_name,strategy,sampling_ratio)
    cat("This is the",i,"bootstrap iteration","\n")
    outcome_sampled = outcome_subsampling_for_binary_trait(outcome,target_outcome_name,sampling_ratio)
    expr = data.frame(expr)
    # Here, after apply the inner_join function, the column name of outcome changed to "X30690.0.0"(original:30690-0.0)
    df_outcome_expr = inner_join(outcome_sampled, expr, by="FID")   #merge(outcome, expr, by="FID")
    #print(  (dim(df_outcome_expr))  )
    #***********************************************************************
    # Revised by: yinly, Date: March 6, 2020
    # Description: preselect a set of genes that are related to the outcome
    #***********************************************************************
    df_outcome_expr_PCs = inner_join(df_outcome_expr,PCs,by ="FID")
    # Commented by yinly, July 9, 2020
    # Revised by yinly, June 15, no need to perform feature selection after removing rankNorm
    df_outcome_expr_PCs_Selected <- df_outcome_expr_PCs
    #df_outcome_expr_PCs_Selected = univariate_gene_selection(df_outcome_expr_PCs,p_thres,outcome_index)
    #print(df_outcome_expr_PCs_Selected)
    cat("The feature dimension for selected expr_outcome_PC dataset is",dim(df_outcome_expr_PCs_Selected),"\n")
    #*****************
    # count the no. of genes available
    #******************
    grep_res = grep("ENSG", colnames(df_outcome_expr_PCs_Selected) )
    start.ind = min(grep_res)
    end.ind = max(grep_res)
    no_genes = end.ind - start.ind + 1
    ##start.ind = index for which we start to have gene expressions
  
    #********************************
    # Residualized X
    #*********************************
    resid.mat = matrix(nrow = nrow(df_outcome_expr_PCs_Selected), ncol= no_genes )
    attach(df_outcome_expr_PCs_Selected)
    for (k in 1:no_genes) {
      lm.obj = lm(df_outcome_expr_PCs_Selected[,start.ind+ k-1] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
      resid.mat[,k] = lm.obj$resid
    }
  
    colnames(resid.mat) <- colnames(df_outcome_expr_PCs_Selected)[start.ind:end.ind]
  
    #********************************
    # Residualized Y
    #*********************************
    # Revised by yinly on Jan 7,2020,the inner_join function changed the column name of the target_outcome_name
    # # Reset the value of
    # outcome_index = "X30690.0.0"
    lm.obj_y = lm(df_outcome_expr_PCs_Selected[,outcome_index] ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
    resid_y = lm.obj_y$resid
    detach(df_outcome_expr_PCs_Selected)
  
    dfnew = data.frame(outcome= resid_y,  resid.mat)
    #save(dfnew, file= "df_new_original_wholeBlood_UKBB.Rdata")
    #*************************
    # Rank-based Inverse normal transformaton
    #************************
    # Revised by:yinly, June 15, 2021
    # Description: no need to perform rankNorm(change distribution of binary traits)
    #dfnew_rankNorm = apply(dfnew, 2, rankNorm)
    dfnew_rankNorm <- dfnew
    library(coop)
    #**************************************
    #  IMPORTANT: the following expression corr matrix can be re-used for other outcomes (given we are working on the SAME tissue)
    #**************************************
    expr_corMat = pcor(dfnew_rankNorm[,-1])
    # save(expr_corMat, file="correl_mat_expr_rankNorm_UKBB_Lung.Rdata")
    #now add back the correlation between the outcome and each feature
    #cf https://stackoverflow.com/questions/20410768/how-to-correlate-one-variable-to-all-other-variables-on-r
    cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
    precompute_corMat = cbind(cor_with_outcome, expr_corMat)
    precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat )
    library(pcalg)
    t1 = proc.time()
    # Revise by yinly:
    # Date: Feb 12,2020
    # Description: set more stringent alpha for binary traits
    pcSimple.fit <- PCSelect_Parallel(y=dfnew_rankNorm[,1],
                                    dm=dfnew_rankNorm[,-1] ,
                                    method = c("parallel"),
                                    mem.efficient = FALSE,
                                    num_workers = 20,
                                    #alpha=0.001,
                                    alpha=0.001,
                                    corMat = precompute_corMat,
                                    max_ord=3,
                                    corMethod = "standard",
                                    verbose = TRUE, directed = TRUE)
    proc.time()-t1
    genes = names(pcSimple.fit$G)
    assoc =  as.numeric(as.vector(unlist(pcSimple.fit$G)))
    zMin = as.vector(unlist(pcSimple.fit$zMin))
    #PCSelect_Result = as.data.frame(cbind(genes,assoc,zMin))
    PCSelect_Result = data.frame(genes = genes,assoc = assoc,zMin = zMin)
    #colnames(PCSelect_Result) = c("genes","assoc","zMin")
    qlambda_sum = qlambda_sum + length(which(PCSelect_Result$assoc==1))
    #PCSelect_Result = PCSelect_Result[order(PCSelect_Result$genes),]
    if(i==1){
      assoc_mat = PCSelect_Result[,1:2]
      zMin_mat = PCSelect_Result[,c(1,3)]
    }else{
      # After ordering the results based on gene names, all bootstrapped results should have the same gene ordering
      # if inner_join function is too time consuming, we can directly use cbind, but need to confirm first
      assoc_mat = merge(assoc_mat, PCSelect_Result[,1:2], by="genes",all=TRUE,sort=FALSE)
      zMin_mat = merge(zMin_mat,PCSelect_Result[,c(1,3)],by="genes",all=TRUE,sort=FALSE)
    }
  }
  colnames(assoc_mat)[2:(iter+1)] = seq(1,iter)
  colnames(zMin_mat)[2:(iter+1)] = seq(1,iter)
  #save(assoc_mat,file = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Bootstrap_Results/Cov19_assoc_mat.Rdata")
  #save(zMin_mat, file = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Bootstrap_Results/Cov19_zMin_mat.Rdata")
  PCSelect_Result_final = matrix(0,nrow=nrow(assoc_mat),ncol=3)
  colnames(PCSelect_Result_final) = c("genes","assoc","zMin")
  PCSelect_Result_final[,"genes"] = as.vector(unlist(assoc_mat$genes))

  # only keep value columns for assoc_mat and zMin mat
  assoc_all = as.matrix(assoc_mat)
  assoc_all = assoc_all[,2:ncol(assoc_all)]
  class(assoc_all) = "numeric"

  zMin_all = as.matrix(zMin_mat)
  zMin_all = zMin_all[,2:ncol(zMin_all)]
  class(zMin_all) = "numeric"
  for(j in 1:nrow(assoc_all)){
    # the cutoff value for the bootstrap strategy is 60
    #if(sum(as.numeric(assoc_mat[j,2:ncol(assoc_mat)]))>iter*sampling_cutoff){
    if(sum(assoc_all[j,],na.rm = TRUE)>iter*sampling_cutoff){ 
      PCSelect_Result_final[j,2] = 1
      assoc_col_index = which(assoc_all[j,]==1)  
    } else{
      # Here we need to confirm whether non-selected genes are considered non-causal, if so, we need to treat them as 0 as well,zMin as 0?
      PCSelect_Result_final[j,2] = 0
      assoc_col_index = which(assoc_all[j,]==0)
    }
    zMin_ave = mean(as.numeric(zMin_all[j,assoc_col_index]))
    PCSelect_Result_final[j,3] = zMin_ave
  }
  # Revised by yinly, May 24, 2021
  # Description: the original calculation for qlambda is wrong, already correct it in this version
  #qlambda = qlambda_sum/nrow(assoc_mat)
  qlambda = qlambda_sum/iter
  cat("qlambda for",target_outcome_name,"is",qlambda,"\n")
  p = nrow(assoc_mat)
  cat("The number of genes(p) for",target_outcome_name,"is",p,"\n")
  EV = (1/(2*sampling_cutoff - 1)) * ((qlambda**2)/p)
  cat("The number of falsely identified causal genes(E(V)) for",target_outcome_name,"is",EV,"\n")
  TypeI_Error = EV/p
  cat("The typeI error for",target_outcome_name,"is",TypeI_Error,"\n")
  FDR = EV/qlambda
  cat("The FDR for",target_outcome_name,"is",FDR,"\n")
  PCSelect_Result_final_list[[1]] = PCSelect_Result_final
  PCSelect_Result_final_list[[2]] = qlambda
  PCSelect_Result_final_list[[3]] = p
  PCSelect_Result_final_list[[4]] = EV
  PCSelect_Result_final_list[[5]] = TypeI_Error
  PCSelect_Result_final_list[[6]] = FDR
  PCSelect_Result_final_list[[7]] = assoc_mat
  PCSelect_Result_final_list[[8]] = zMin_mat
  return(PCSelect_Result_final_list)
}

setwd("/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/Order3_Results/")
target_outcome_name = "CAD"
outcome_index = "CAD"
outcome_name = "CAD"
#strategy = "First"
# Use new extracted phenotypes
#var_thres = 0.01
# Revised by yinly, May 22, 2021
# Description: update the imputed gene expression file and extracted phenotype file accordingly
expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Whole_Blood_CIS_Trans_expression.txt"
pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_30April.txt"
#pheno_dir = "/exeh_3/yinly/UKBB/01_Data/UKBB_Cov19_with_FID.txt"
#___  END OF INPUT__________________________________________
  
expr = fread(expression_file_dir )
#************************************
# screen away genes with very low variance (this will lead to NA entries after scaling in huge
#***********************************
# Commented by yinly, June 25, 2021
# Description: no need to filter expression data based on variances, since the imported data is already filtered!
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
#**********************************************************************
# Revised by: yinly
# Date: July 9, 2020
# Description: select cases and controls based on different strategies
# Cov19: 0(controls),1 non-hospitalized cov19 patients,2 hospitalized cov19 patients
# 1) Hopspitalized cov19 patients vs non-hospitalized cov19 patients
# 2) All cov19 patients vs all controls
# 3) Hospitalized patients vs all other subjects
#*********************************************************************
#outcome = get_outcome_samples_for_cov19(outcome,outcome_index,strategy)

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
colnames(PCs)[1] <- "FID"
colnames(PCs)[2:11] <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")

sampling_ratio = 0.6
# To increase computational efficiency, we set the iter to 50
iter = 50
sampling_cutoff = 0.6
p_thres = 0.2
pcSimple.fit = Get_Bootstrapped_Gene_Outcome_Network(expr,outcome,target_outcome_name,sampling_ratio,iter,PCs,sampling_cutoff,p_thres)
pcSimple_prefix = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Stability_Selection_Results/"
#pcSimple_name = paste0("UKBB_",outcome_name,"_Whole_Blood_",strategy,".Rdata")
pcSimple_name = paste0("UKBB_",outcome_name,"_Whole_Blood.Rdata")
save(pcSimple.fit,file=paste0(pcSimple_prefix,pcSimple_name))
#fwrite(as.data.frame(pcSimple),paste0(pcSimple_prefix,pcSimple_name),sep="\t")
