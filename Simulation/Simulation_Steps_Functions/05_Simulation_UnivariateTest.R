#in this function we can have the correlation of row data and imputed data.
Univariate_Test<-function(gene_index,samples_data,samples_data_iptg)
{
  outcome = samples_data[,dim(samples_data)[2]] # phenotype or outcome is the last column in this sample matrix
  
  
  univariate_test_Matrix = matrix(nrow = length(gene_index), ncol = 5)
  
  count = 1
  for(gene in gene_index)
  {
    print(gene)
    gene_rowdata = samples_data[,gene]
    gene_imputedata = samples_data_iptg[,gene]
    
    lm_gene_rowdata = lm(outcome~gene_rowdata)
    lm_gene_imputedata = lm(outcome~gene_imputedata)
    
    lm_gene_rowdata_pval = summary(lm_gene_rowdata)$coefficients[2,4]
    lm_gene_imputedata_pval = summary(lm_gene_imputedata)$coefficients[2,4]
    
    # univariate_test_Matrix[count,1] = gene # write gene_index into this column 
    univariate_test_Matrix[count,1] = count
    univariate_test_Matrix[count,2] = unname(coefficients(lm_gene_rowdata)[2])
    univariate_test_Matrix[count,3] = unname(coefficients(lm_gene_imputedata)[2])
    univariate_test_Matrix[count,4] = lm_gene_rowdata_pval
    univariate_test_Matrix[count,5] = lm_gene_imputedata_pval
    
    count = count + 1
  }
  return(univariate_test_Matrix)
}


