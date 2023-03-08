library(data.table)
library(dplyr)

get_univariate_test_results <-function(trait,tissue,target_outcome_name,outcome_index){
    expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Whole_Blood_CIS_Trans_expression.txt"
    predicted_pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/02_Multiple_Traits/Genetically_Predicted_HbA1c.csv"
    #pheno_dir = "/exeh_3/yinly/UKBB/01_Data/UKBB_Traits_with_ID_July29.csv"
    #pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_with_FID_updated.csv"
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

    #************************************
    # Load genetically predicted pheno and combine with gene expression data
    #************************************
    predicted_pheno = fread(predicted_pheno_dir)
    colnames(predicted_pheno) = c("FID","ENSG_HbA1c")
    expr_new = inner_join(expr, predicted_pheno,by="FID")

    
  
    #********************************
    # Extract outcome trait - Just take BMI as an example; code can be searched from the HTML file
    #********************************
    pheno = fread(pheno_dir)
    outcome = cbind(pheno[,1], pheno[,..target_outcome_name])
    outcome  = data.frame(   na.omit(outcome)   )

    expr = data.frame(expr_new)
    # Here, after apply the inner_join function, the column name of outcome changed to "X30690.0.0"(original:30690-0.0)
    df = inner_join(outcome, expr, by="FID")   #merge(outcome, expr, by="FID")
    print(  (dim(df))  )
  
    pheno = data.frame(pheno)
    PC_start_ind = which( colnames(pheno) == "X22009.0.1")
    PCs = cbind(pheno[,"FID"], pheno[,PC_start_ind:(PC_start_ind+10-1)]  )
    
    colnames(PCs)[1] <- "FID"
    colnames(PCs)[2:11] <- c("PC1", "PC2", "PC3","PC4","PC5","PC6","PC7","PC8","PC9", "PC10")
    df = inner_join(df,PCs,by ="FID")
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
    for (i in 1:no_genes) {
      # Commented by:yinly
      # Perform univariate test by linear regression to select significant genes
      lm.obj = lm(df[,outcome_index] ~ df[,start.ind+ i-1] + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 )
      #resid.mat[,i] = lm.obj$resid
      univariate.mat[i,1] = colnames(df)[start.ind + i - 1]
      univariate.mat[i,2] = summary(lm.obj)$coeff[2,1]
      univariate.mat[i,3] = summary(lm.obj)$coeff[2,4]
    }
  
    colnames(univariate.mat) = c("Genes","Estimates","Pvalues")
    #tissue = "Whole_Blood"
    Univariate_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Univariate_Results_Update/"
    Univariate_filepath = paste0(Univariate_prefix,tissue,"_HbA1c_",trait,"_Univariate_Linear.csv")
    fwrite(as.data.frame(univariate.mat),Univariate_filepath)
}

Dict_filepath = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Whole_Blood_associated_traits_update.csv"
Dict = as.matrix(fread(Dict_filepath))
#index <- c(1,4)
#for(j in 1:2){
  #i <- index[j]
  #i = 21
  #trait = "Cov19"
  #tissue = "Whole_Blood"
  #target_outcome_name = "Cov19"
  #outcome_index = "Cov19"
  i <- 4
  trait = Dict[i,1]
  tissue = Dict[i,2]
  target_outcome_name = Dict[i,3]
  outcome_index = Dict[i,4]
  cat("Univariate test for the",i,"trait",trait,"with",tissue,"as causal tissue","\n")
  get_univariate_test_results(trait,tissue,target_outcome_name,outcome_index)
  cat("Univariate test for",trait,"with",tissue,"as causal tissue is completed!","\n")
#}
