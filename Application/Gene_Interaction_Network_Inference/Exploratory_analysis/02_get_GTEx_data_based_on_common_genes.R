#*************************************************************************************
# Written by:yinly, May 21, 2020
# Description: include more genes for the GTEx network(Causal and noncausal genes)
#*************************************************************************************

#*************************************
#  load libraries
#***************************************
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(dplyr)
library(pcalg)
library(coop)
library(parallel)

Get_Residual_mat_for_GTEx_Network <- function(trait,tissue_type,q_thres){
  #***************************************************************************************
  # Step 1: load tissue-specific residualized gene expression 
  #***************************************************************************************
  tissue_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Results/" # directory for the adjusted GTEx tissue-specific expression data
  tissue_resid_name = "Resid_Whole_Blood_GTEx.Rdata" # filename for the adjusted tissue-specific expression data
  load(paste0(tissue_file_prefix,tissue_resid_name))
  
  #***************************************************************************************
  # Step 2: load tissue traits results, select the top 100 genes 
  # if no. associated genes > 100, we select the top 100 with largest zMin
  #***************************************************************************************
  tissue_trait_file_prefix = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/Order3_Results_Update/" # directory for the inferred gene-outcome causal networks for studied traits
  #trait = "CAD"
  #tissue_type = "Adipose_Visceral_Omentum"
  #if(tissue_type =="Adipose_Visceral_Omentum"){
  #  tissue_trait_file_name = paste0(tissue_trait_file_prefix,"UKBB_",trait,".Rdata")
  #}else{
  #  tissue_trait_file_name = paste0(tissue_trait_file_prefix,"UKBB_",trait,"_",tissue_type,".Rdata")
  #}
  tissue_trait_file_name = paste0(tissue_trait_file_prefix,"UKBB_",trait,"_",tissue_type,".Rdata")
  load(tissue_trait_file_name)
  # Get selected genes from pcSimple.fit
  genes_names = names(pcSimple.fit$G)
  assoc =  as.vector(unlist(pcSimple.fit$G))
  zMin = as.vector(unlist(pcSimple.fit$zMin))
  tissue_trait_result = data.frame(genes_names,assoc,zMin,stringsAsFactors=FALSE)
  
  #*************************************************************************************
  # Revised by yinly, June 9,2021
  # Description: no need to order PCSelect results based on zMin. If so, we only need to 
  # infer one gene-gene network for each tissue. This would be more efficient. Previously 
  # We ordered this because there is inconsisency between binary traits. We have solve this
  # problem by remove rankNorm, thus, we don't have to reorder this anymore.
  #****************************************************************************************
  #UKBB_selected = tissue_trait_result[tissue_trait_result$assoc==TRUE,]
  #UKBB_selected_sort = tissue_trait_result[order(-tissue_trait_result$zMin),]
  UKBB_selected_sort <- tissue_trait_result
  #*************************************************************************************
  # Revised by:yinly, May 13, 2020
  # Description: to get the top 100 genes based on zMin, since there are some discrepancy 
  # between genes in GTEx and UKBB Imputed expression, we select the top 100 among shared
  # genes if no.genes > 100
  #****************************************************************************************
  # Commented by yinly, May 21,2020
  # If we change the strategy for gene selection, we should change the following few lines
  #****************************************************************************************
  Genes_names_all = colnames(resid_mat) # Get GTEx gene names
  Genes_names_UKBB = tissue_trait_result[,1] # Get UKBB gene names from PCSelect
  Genes_intersect = intersect(Genes_names_UKBB,Genes_names_all) # Get shared genes between GTEx and UKBB imputed tissues
  Genes_shared = tissue_trait_result[(tissue_trait_result[,1] %in% Genes_intersect),] # Remove nonshared genes from UKBB
  # Revised by yinly, June 9,2021
  # Description: no need to order genes based on zMin
  #Genes_shared_sort = Genes_shared[order(-Genes_shared$zMin),] # Order shared genes based on zMin in descending order
  Genes_shared_sort <- Genes_shared

  Genes_selected = Genes_shared_sort$genes_names
  
  Genes_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_All_Selected_Genes_CIS_Trans_Update/" # directory for the files recording all selected genes for each studied tissue-specific trait pairs
  GTEx_selected_genes = Genes_shared_sort[(Genes_shared_sort[,1]%in%Genes_selected),]
  cat("The dimension of selected genes for",tissue_type,"related trait",trait,"is",dim(GTEx_selected_genes),"\n")
  fwrite(GTEx_selected_genes,paste0(Genes_file_prefix,tissue_type,"_",trait,"_Selected_Genes.csv"))
  resid_mat_selected = resid_mat[,Genes_selected]
  resid_mat_selected_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/Residmat_All_Update/" # directory used to store the adjusted GTEx tissue-specific expression profile for common genes between GTEx and imputed tissue-specific expression profiles for subjects in UK Biobank
  resid_mat_selected_filepath = paste0(resid_mat_selected_prefix,tissue_type,"_",trait,"_resid_mat_selected.Rdata")
  save(resid_mat_selected,file=resid_mat_selected_filepath)
}

setwd("/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/")
#********************************************************************
# The following 5 lines of codes are for testing only
#********************************************************************
# trait = "CAD"
# tissue_type ="Prostate"
# cat("This is the",tissue_type,"tissue for",trait,"\n")
# Get_Residual_mat_for_GTEx_Network(trait,tissue_type)
# cat("The",tissue_type,"specific GTEx network for",trait,"is completed!","\n")
#tissue_type = "Adipose_Visceral_Omentum"
traits_dict = as.matrix(fread("Whole_Blood_associated_traits_updated.csv")) # csv file records whole blood associated traits

for(i in 1:nrow(traits_dict)){
  trait = traits_dict[i,1]
  tissue_type = traits_dict[i,2]
  cat("This is the",tissue_type,"tissue for",trait,"\n")
  Get_Residual_mat_for_GTEx_Network(trait,tissue_type,0.1)
  cat("The",tissue_type,"specific GTEx network for",trait,"is completed!","\n")
}
#trait <- "Breast_cancer"
#tissue_type <- "Breast_Mammary"
#cat("This is the",tissue_type,"tissue for",trait,"\n")
#Get_Residual_mat_for_GTEx_Network(trait,tissue_type,0.1)
#cat("The",tissue_type,"specific GTEx network for",trait,"is completed!","\n")
