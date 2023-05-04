
#***********************************************************
#  The aim of this code is to prepare the gene expression matrix, controlling for possible batch effects (and other confounding factors) 
#***********************************************************
#************************
#   INPUT 
#************************

TPM_thres = 0.2  ##can change; genes with median TPM lower than this value are removed 
low_var_thres = 0.1  ## genes with variance in expression lower than this value will be removed
tissue_type_to_study = "Whole Blood"


#*************************************
#  pc algorithm applied to GTEx 
#***************************************
setwd("/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Results/")
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(dplyr)
library(pcalg)
library(coop)
library(parallel)



load("/exeh_4/sohc/network_GWAS/GTEx_RNA_seq_transposed_df.Rdata") # transposed GTEx RNA-seq data downloaded from GTEx(https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz)
attribute = fread("/exeh_4/sohc/network_GWAS/GTEx_v7_Annotations_SampleAttributesDS.txt", header=TRUE) # downloaded data from GTEx(https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt)
##filter away genes with low TPM
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

##WRONG
##attribute$SUBJID <- substr(attribute$SAMPID, 1,10) 
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
subject_pheno = fread("/exeh_4/sohc/network_GWAS/GTEx_v7_Annotations_SubjectPhenotypesDS.txt") # downloaded phenotype file form GTEx(https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt)
SUBJID_dot <- NULL
for (i in 1:nrow(subject_pheno))  { 
   SUBJID_dot[i] <- gsub("-", '.', subject_pheno[,"SUBJID"][i])  
}
subject_pheno$SUBJID <- SUBJID_dot

finaldf3 = merge(finaldf3,  subject_pheno, by = "SUBJID")

# library( DMwR)
# finaldf3_imp = knnImputation(data= finaldf3, k = 10, scale = T, meth = "median") 

#******************************************
# PCA analysis to remove batch effects (or other hidden confounders) in gene co-expression 
# https://www.biorxiv.org/content/biorxiv/early/2018/10/05/202903.full.pdf
#*******************************************************************
library(gmodels)
#pca_res = fast.prcomp(finaldf3[,3:(ind_start_attribute-1)] ) 
#save(pca_res, file= "/exeh_4/sohc/network_GWAS/R_code_BayesianNet/PCA_result_GTEx_expression.Rdata") 
load("/exeh_4/sohc/network_GWAS/R_code_BayesianNet/PCA_result_GTEx_expression.Rdata") # PCA analysis results for the GTEx expression.

PCA = data.frame(pca_res$x[,1:10]) ##Cf https://rpubs.com/skydome20/R-Note7-PCA
finaldf3 <- data.frame(finaldf3, PCA) 




#******************************
# The expression of each gene is regressed onto tissue type
# This is because two genes may appear correlated merely because they are both over-/down-expressed in a tissue
# in fact one gene does not CAUSE expression change in another
#  ie tissue type can be a confounding factor 
# ALTERNATIVELY, you can specify one tissue 
#***********************************
tissue_type = finaldf3$SMTSD
print(unique(tissue_type))
#*****************************************
# (IMPORATNT)
# IF need to retain only one tissue type, change this line
#******************************************
final_df_SingleTissue = finaldf3[finaldf3$SMTSD==tissue_type_to_study,]  

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
save(resid_mat, file= "Resid_Whole_Blood_GTEx.Rdata") 
