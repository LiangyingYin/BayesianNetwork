#******************************************************************************************
# Written by yinly, June 24, 2021
# Description: after removing rankNorm for the binary traits for PCSimple algorithm,
# there are only few genes identified as directly causal to our studied phenotypes under
# the current cutoff, we could change the strategy for gene selection for residmat, i.e.,
# if the no. of directly causal genes is less than 20, then we select the top 20 genes with 
# the largest zMins, this also applied for the gene selection strategy for univariate test
#*******************************************************************************************
library(data.table)
#**********************************************************************
# para@trait: studied trait
# para@tissue_type: studied tissue 
# para@p_thres: p-value cutoff used for univariate gene selection
#**********************************************************************
Get_Residual_mat_for_GTEx_Network <- function(trait,tissue_type,p_thres){
  #***************************************************************************************
  # Step 1: load tissue-specific residualized gene expression 
  #***************************************************************************************
  tissue_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Results/"
  tissue_resid_name = "Resid_Lung_GTEx.Rdata"
  load(paste0(tissue_file_prefix,tissue_resid_name))
  
  #***************************************************************************************
  # Step 2: load tissue traits results, select the top 100 genes 
  # if no. associated genes > 100, we select the top 100 with largest zMin
  #***************************************************************************************
  tissue_trait_file_prefix = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/Order3_Results_Update/"
  tissue_trait_file_name = paste0(tissue_trait_file_prefix,"UKBB_",trait,"_",tissue_type,".Rdata")
  load(tissue_trait_file_name)
  # Get selected genes from pcSimple.fit
  genes_names = names(pcSimple.fit$G)
  assoc =  as.vector(unlist(pcSimple.fit$G))
  zMin = as.vector(unlist(pcSimple.fit$zMin))
  tissue_trait_result = data.frame(genes_names,assoc,zMin,stringsAsFactors=FALSE)
  UKBB_sort <- tissue_trait_result[order(-tissue_trait_result$zMin),]
  UKBB_selected = tissue_trait_result[tissue_trait_result$assoc==TRUE,]
  UKBB_selected_sort = UKBB_selected[order(-UKBB_selected$zMin),]
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
  Genes_intersect = intersect(Genes_names_UKBB,Genes_names_all) # Get shared genes between GTEx and UKBB imputed tisues
  Genes_shared = tissue_trait_result[(tissue_trait_result[,1] %in% Genes_intersect),] # Remove nonshared genes from UKBB
  Genes_shared_sort = Genes_shared[order(-Genes_shared$zMin),] # Order shared genes based on zMin in descending order
  
  # Load significant genes based on univariate test
  Univariate_prfix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Univariate_Results_All_Rerun/"
  Univariate_filepath = paste0(Univariate_prfix,tissue_type,"_",trait,"_Univariate_Linear.csv")
  Univariate = as.matrix(fread(Univariate_filepath))
  Qvalues = p.adjust(as.numeric(Univariate[,"Pvalues"]),method="BH")
  Univariate = cbind(Univariate,Qvalues)
  # Revised by yinly, June 24, 2021
  # Description: only keep genes shared with GTEx for further analysis
  Univariate <- Univariate[Univariate[,1] %in% Genes_shared_sort[,1],]
  Univariate_sort = Univariate[order(as.numeric(Univariate[,"Pvalues"])),]
  Univariate_sig = Univariate_sort[(as.numeric(Univariate_sort[,"Pvalues"])<p_thres),]
  
  # Get shared important genes derived from PCSelect(Selected genes in UKBB imputed tissues may not exist in GTEx)
  UKBB_selected_genes_shared = Genes_shared_sort[Genes_shared_sort$assoc==TRUE,1]
  
  # General strategy for gene selection.
  # if the no. selected genes for UKBB<100,select same number of genes based on univariate test.
  # We need to consider the special cases where only one gene remaining for the univariate test
  # Revised by yinly, June 24, 2021
  # Description: if the no. of directly causal genes is less than 20, then we selected the top 20
  # with the largest zMin, this also apply for the univariate test
  # Revised by yinly: August 11, 2021
  # Description: we select the top 150 directly causal genes for the analysi
  if(length(UKBB_selected_genes_shared)<20){
      # no. directly causal genes <20, then selecte the top 20 with largest zMins and top 20 with smallest p for univariate test
      UKBB_selected_genes <- Genes_shared_sort[1:20,1]
      #Univariate_genes <- Univariate_sort[1:20,1]
      cat("The number of directly causal genes is less than 20, we select top 20 with largest zMins and top 20 with smallest pvalue for univariate test! \n")
  } else if(length(UKBB_selected_genes_shared)<150){
      UKBB_selected_genes <- UKBB_selected_genes_shared
      #Univariate_genes <- Univariate_sort[1:length(UKBB_selected_genes_shared),1]
      cat ("The number of directly causal genes is between 20 and 150, all of them plus the same no. of significant genes from univariate test are selected for analysis! \n")
  } else{
      UKBB_selected_genes <- UKBB_selected_genes_shared[1:150]
      #Univariate_genes <- Univariate_sort[1:100,1]
      cat ("The number of directly causal genes is more than 150, the top 150 with largest zMins plus the same no. of significant genes from univariate test are selected for analysis! \n")
  }
  #Genes_selected = union(UKBB_selected_genes,Univariate_genes)
  # Revised by yinly, August 11,2021
  # Only include the genes from PCSimple
  Genes_selected <- UKBB_selected_genes
  
  # Save selected important genes derived from the GTEx and UKBB imputed tissue respectively.
  Genes_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_Selected_Genes_Update/"
  GTEx_selected_genes = Genes_shared_sort[(Genes_shared_sort[,1]%in%Genes_selected),]
  cat("The dimension of selected genes for",tissue_type,"related trait",trait,"is",dim(GTEx_selected_genes),"\n")
  fwrite(GTEx_selected_genes,paste0(Genes_file_prefix,tissue_type,"_",trait,"_Selected_Genes.csv"))
  resid_mat_selected = resid_mat[,Genes_selected]
  resid_mat_selected_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Result_Selected_Update/"
  resid_mat_selected_filepath = paste0(resid_mat_selected_prefix,tissue_type,"_",trait,"_resid_mat_selected.Rdata")
  save(resid_mat_selected,file=resid_mat_selected_filepath)
}


setwd("/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/")
#********************************************************************
# The following 5 lines of codes are for testing only
#********************************************************************
# trait = "Breast_cancer"
# tissue_type ="Breast_Mammary"
# cat("This is the",tissue_type,"tissue for",trait,"\n")
# Get_Residual_mat_for_GTEx_Network(trait,tissue_type,0.001)
# cat("The",tissue_type,"specific GTEx network for",trait,"is completed!","\n")

traits_dict = as.matrix(fread("Lung_associated_traits_updated.csv"))
#for(i in 1:nrow(traits_dict)){
for(i in 1:2){
 trait = traits_dict[i,1]
 tissue_type = traits_dict[i,2]
 cat("This is the",tissue_type,"tissue for",trait,"\n")
 Get_Residual_mat_for_GTEx_Network(trait,tissue_type,0.001)
 cat("The",tissue_type,"specific GTEx network for",trait,"is completed!","\n")
}
