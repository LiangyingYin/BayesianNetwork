
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(dplyr)
library(pcalg)
library(coop)
library(parallel)
library(glasso)
library(qgraph)
library(gmodels)
library(glassoFast)
library(Matrix)

#*************************************************************************
#  The aim of this function is to prepare the gene expression matrix, 
#  controlling for possible batch effects (and other confounding factors) 
#*************************************************************************

get_Residualized_GTEx <- function(TPM_thres,low_var_thres,tissue_type_to_study,resid_file_prefix){
  load("/exeh_4/sohc/network_GWAS/GTEx_RNA_seq_transposed_df.Rdata") 
  attribute = fread("/exeh_4/sohc/network_GWAS/GTEx_v7_Annotations_SampleAttributesDS.txt", header=TRUE)
  #****************************************************
  #  here filtering genes with low overall expression may be justified as such genes are unlikely to be biologically important 
  #  ie they are unlikely to be important confounders and unlikey to be causal genes for diseases
  #***************************************************** 
  median_TPM = apply(finaldf[,-1], 2, median) 
  mean_TPM = colMeans(finaldf[,-1]) 
  
  ind_highTPM = which(median_TPM>= TPM_thres)
  finaldf = finaldf[, c(1,(ind_highTPM+1)) ]
  #********************************************
  #     log transformed the TPM metrics
  #  cf: https://support.bioconductor.org/p/104943/
  #********************************************
  SAMPID = finaldf[,"variable"]
  finaldf2 = data.frame(SAMPID,  log2(finaldf[,-1]+1) ) 
  
  #***************************************
  #  filter genes with low variance in expression 
  # as two genes with low varinace will be highly correlated cf: https://bioconductor.org/packages/release/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html
  #*****************************************
  var_gtex = apply(finaldf2[,-1], 2, var)
  ind_high_var = which(var_gtex >= low_var_thres )  
  finaldf2 = finaldf2[, c(1,(ind_high_var+1)) ]
  
  SAMPID_dot <- NULL
  #replace - with . to make the SAMPID identical
  for (i in 1:nrow(attribute))  { 
    SAMPID_dot[i] <- gsub("-", '.', attribute[,"SAMPID"][i])  
  }
  attribute$SAMPID <- SAMPID_dot
  
  SUBJID = NULL
  splitted  <- strsplit(attribute$SAMPID, "[.]")
  for (i in 1:length(attribute$SAMPID) ) {
    SUBJID[i] = paste( splitted[[i]][1] , splitted[[i]][2], sep=".")
  }
  
  attribute$SUBJID  <- SUBJID
  finaldf3 = merge(finaldf2, attribute, by="SAMPID")
  
  #**********************************
  # read in the subject phenotype file
  #********************************
  subject_pheno = fread("/exeh_4/sohc/network_GWAS/GTEx_v7_Annotations_SubjectPhenotypesDS.txt")
  SUBJID_dot <- NULL
  for (i in 1:nrow(subject_pheno))  { 
    SUBJID_dot[i] <- gsub("-", '.', subject_pheno[,"SUBJID"][i])  
  }
  subject_pheno$SUBJID <- SUBJID_dot
  finaldf3 = merge(finaldf3,  subject_pheno, by = "SUBJID")
  
  #******************************************
  # PCA analysis to remove batch effects (or other hidden confounders) in gene co-expression 
  # https://www.biorxiv.org/content/biorxiv/early/2018/10/05/202903.full.pdf
  #*******************************************************************
  
  load("/exeh_4/sohc/network_GWAS/R_code_BayesianNet/PCA_result_GTEx_expression.Rdata")
  
  PCA = data.frame(pca_res$x[,1:10]) ##Cf https://rpubs.com/skydome20/R-Note7-PCA
  finaldf3 <- data.frame(finaldf3, PCA) 
  #******************************
  # The expression of each gene is regressed onto tissue type
  # This is because two genes may appear correlated merely because they are both over-/down-expressed in a tissue
  # in fact one gene does not CAUSE expression change in another
  # ie tissue type can be a confounding factor 
  # ALTERNATIVELY, you can specify one tissue 
  #***********************************
  tissue_type = finaldf3$SMTSD
  #print(unique(tissue_type))
  #*****************************************
  # (IMPORATNT)
  # IF need to retain only one tissue type, change this line
  #******************************************
  final_df_SingleTissue = finaldf3[finaldf3$SMTSD==tissue_type_to_study,]  
  print(paste0(tissue_type_to_study,":",nrow(final_df_SingleTissue)," samples"))
  print(dim(final_df_SingleTissue))
  print(length(which(is.na(final_df_SingleTissue))))
  if(length(which(is.na(final_df_SingleTissue)))>0){
     print(which(is.na(final_df_SingleTissue))/ncol(final_df_SingleTissue))
  }
  ind_start_attribute = which(colnames(final_df_SingleTissue)=="SMATSSCR")
  no_genes = ind_start_attribute-3
  
  #resid_mat = matrix(ncol= no_genes, nrow= nrow(final_df_SingleTissue)-103 )  ##missing DTHHRDY 103 ##need to change this no. if add other covariates 
  resid_mat = matrix(ncol= no_genes, nrow= nrow(final_df_SingleTissue) ) 
  
  attach(final_df_SingleTissue) 
  for (i in 1:no_genes) {
    lm_obj = lm( final_df_SingleTissue[,i+2]~  PC1 + PC2 +PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                   #eigenvec$PC1 + eigenvec$PC2 + eigenvec$PC3 + eigenvec$PC4 + eigenvec$PC5 + eigenvec$PC6 + eigenvec$PC7 + eigenvec$PC8 + eigenvec$PC9 + eigenvec$PC10 +
                   #factor(final_df_SingleTissue$SMTSD) + 
                   factor(SEX) + 
                   factor(AGE) + 
                   #factor(DTHHRDY) +  #missing =103
                   #SMATSSCR + ##missing 1245
                   factor(SMGEBTCH) +   #no missing 
                   SMRIN    #no missing
                 #SMTSISCH + ##NA = 952
                 #SMTSPAX  #missing 1245
    )
    resid_mat[,i] = lm_obj$residuals
  }
  detach(final_df_SingleTissue) 
  
  colnames(resid_mat) <- colnames(final_df_SingleTissue)[3:(ind_start_attribute-1)]
  tissue_name = gsub(" - ","_",tissue_type_to_study)
  tissue_name = gsub(" ","_",tissue_name)
  outfile_path = paste0(resid_file_prefix,"Resid_",tissue_name,"_GTEx.Rdata")
  save(resid_mat, file= outfile_path) 
  return(resid_mat)
}


#***************************************************************************************************
# Aim of this function is to employ Graphical lasso (via glassoFast) to learn an independence graph 
# of the gene expressions by applying glasso, nodes that are connected means they are dependent 
# when conditioned on ALL other variables in the graph
#***************************************************************************************************
get_Glasso_Graph <- function(resid_mat,n.cores,tissue_type_to_study,Glasso_file_prefix){
  #_________________________________
  cov_mat = coop::covar(resid_mat)
  t1=proc.time()
  rho_theoretical =  sqrt( log(ncol(resid_mat)) / nrow(resid_mat) )     ## https://arxiv.org/pdf/1403.6752.pdf p.1219 / SILGGM publication https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006369
  print(rho_theoretical)
  glasso_obj = glassoFast(S=cov_mat, 
                          rho= rho_theoretical) 
  
  proc.time()-t1
  tissue_name = gsub(" - ","_",tissue_type_to_study)
  tissue_name = gsub(" ","_",tissue_name)
  outfile_path = paste0(Glasso_file_prefix,"Glasso_obj_",tissue_name,"_GTEx.Rdata")
  save(glasso_obj, file = outfile_path) 
  return(glasso_obj)
}

#***************************************************************************************************
# Computes the Bayesian network for GTEx (tissue-specific), using Glasso output as pre-screening step 
 # Edges that are absent in glasso (ie partial correlation=0 after controlling for all other genes) 
# are also set to be absent in pcalg fit
#***************************************************************************************************
get_Bayesian_Network <- function(alpha,numCores,m.max,resid_mat,glasso_obj,PC_file_prefix){
  #*************************************
  #  pc algorithm applied to GTEx 
  #*************************************
  # pcor is a function from "coop" for fast computation of correlation matrix
  corr_mat = pcor(resid_mat) 
  suffStat <- list(C=corr_mat ,n=nrow(resid_mat)) 
  
  #******************************************************************************
  # enforce symmetricity in the partial correlation matrix obtained from glasso
  #******************************************************************************
  # https://stackoverflow.com/questions/18165320/creating-a-symmetric-matrix-in-r
  FixedGaps_mat = forceSymmetric(glasso_obj$wi)
  FixedGaps_mat[FixedGaps_mat!=0] <- 999 
  FixedGaps_mat[FixedGaps_mat==0] <- 1
  FixedGaps_mat[FixedGaps_mat==999] <- 0
  FixedGaps_mat = matrix( as.logical(FixedGaps_mat), nrow = nrow(FixedGaps_mat) ) 
  
  t1=proc.time()
  pc_obj = pc(suffStat, 
              indepTest=gaussCItest, 
              p=ncol(corr_mat), 
              skel.method="stable.fast",
              fixedGaps = FixedGaps_mat,
              alpha = alpha, 
              numCores = numCores ,
              m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209  
              verbose= FALSE)
  proc.time()-t1
  tissue_name = gsub(" - ","_",tissue_type_to_study)
  tissue_name = gsub(" ","_",tissue_name)
  outfile_path = paste0(PC_file_prefix,"PC_all_genes",tissue_name,"_GTEx.Rdata")
  save(pc_obj, file = outfile_path) 
}
