
Imupte_GeneExpressionData<-function(samples_data,cross_validation_data)
{
  # allocate snps data and gene data randomly
  snps_data = cross_validation_data[,snp_index]  ##the 800 subjects (GTEx)
  gene_data = cross_validation_data[,gene_index] ##the 800 subjects (GTEx)
  snps_data2 = samples_data[,snp_index] #use for prediction
  
  
  samples_data_iptg  = samples_for_impute
  
 
  for (k in gene_index)
  {
    
    gene_k = cross_validation_data[,k] # obtain the gene data in turn 
    
    cvfit =cv.glmnet(snps_data,gene_k,nfolds = cv.glmnet_nfolds)
    
    coef(cvfit, s = "lambda.min")
    
    
    #choosing the best lambda according to cross-validation by using  cv.glmnet
    pred_genes = predict(cvfit, newx = snps_data2[1:sample_num_predict_imputedgene,], s = "lambda.min")
    if (var(pred_genes)==0) {
      
      cvm_min_nonzero = min( cvfit$cvm[cvfit$nzero>0] )  #extract the best cv error at which there is at least one non-zero coeff. 
      lamb = cvfit$lambda[cvfit$cvm==cvm_min_nonzero]
      pred_genes = predict(cvfit, newx = snps_data2[1:sample_num_predict_imputedgene,], s = lamb)
    }
    
  
    samples_data_iptg[,k] = pred_genes[,1:1]  #use the 1st column data temporarily
    
  
  }
  
  return(samples_data_iptg)
  
}

Caculate_Cor_ImputeGeneAndRawGene<-function(gene_index,samples_data,samples_data_iptg,sample_num)
{
  dt_genedata = data.frame()
  row_count = 1
  
  for (g_count in gene_index) # iterate over each gene according to gene_index.
  {

    v1 = samples_data[,g_count]
    
    v2 = samples_data_iptg[,g_count]
    
    # results in comparsion
    MSE_value = MSE(v1,v2)/sample_num 
    cor_value = cor(v1,v2)  # in some cases, cor value is NA, because the standard deviation is zero
    cosine_value = cosine(v1,v2) 
    
    dt_genedata[row_count,1] = paste0("G",g_count)
    dt_genedata[row_count,2] = MSE_value
    dt_genedata[row_count,3] = cor_value
    dt_genedata[row_count,4] = cosine_value
    colnames(dt_genedata)[1] = "Gene_num"
    colnames(dt_genedata)[2] = "MSE_value"
    colnames(dt_genedata)[3] = "cor_value"
    colnames(dt_genedata)[4] = "cosine_value"
    row_count = row_count+1
  }
  
  return(dt_genedata) #print gene imputed data and simulated data
  
}

GeneToPhenotype_Network<-function(samples_data_iptg)
{
  
  library(dplyr) 
  library(data.table)
  
  
  var_thres = 0.01
  target_outcome_name = "X20016.0.0"  #Diastolic blood pressure, automated reading
  
  pheno_dir = "/exeh_4/yaning_feng/04_Simulation/Data/array_phenos.txt"
 
  column3 = s+g+1
  column4 = s+g+p
  
  pheno_data = as.matrix(samples_data_iptg[1:sample_num_predict_imputedgene, column3:column4])
  
 
  gene_data = samples_data_iptg[1:sample_num_predict_imputedgene,gene_index] # extract the gene data according to gene_index
  
  
  dfnew_rankNorm = cbind(pheno_data,gene_data) #鐠嬪啯宕瞘ene閸滃henotype閻ㄥ嫭鏆熼幑顕嗙礉閹跺henotype閺佺増宓侀弨鎯ф躬閸撳秹�??
  colnames(dfnew_rankNorm)[1] <-"outcome"
  
  
  library(coop) 
  
  expr_corMat = pcor(dfnew_rankNorm[,-1])
  
  
  cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
  precompute_corMat = cbind(cor_with_outcome, expr_corMat)
  
  precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat) 
  
  
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
  
  
  pcSimple.fit <- pcSelect_stableVer(y=dfnew_rankNorm[,1], 
                                     dm=dfnew_rankNorm[,-1] , 
                                     corMat = precompute_corMat, 
                                     alpha=0.001, 
                                     max_ord=3, 
                                     corMethod = "standard",
                                     verbose = TRUE, directed = TRUE)  
  
  proc.time()-t1
  file_path = "/exeh/exe4/yaning_feng/BayesianNetwork/04_Simulation/Results/fyntest_simulate.Rdata"
  # save(pcSimple.fit, file=file_path)
  return(pcSimple.fit)
}


