
library(data.table)
setwd("/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/")
source("src/GTEx_BN_Inference.R")

tissue_list = as.matrix(fread("list_of_tissues_GTEx.csv"),header = TRUE)
tissues = tissue_list[,2]
#********************************************************************************************************
# Define parameter values for BN inference
#********************************************************************************************************
# Define parameter values for the 1st step:get_Residualized_GTEx
TPM_thres = 0.2
low_var_thres = 0.1
resid_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Results/"

# Define parameter values for the 2nd step:get_Glasso_Graph
n.cores = 15
Glasso_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Glasso_Results/"

# Define parameter values for the 3rd step:get_Bayesian_Network
alpha = 0.01   #alpha (p-value) threshold for 
numCores = 16  #number of cores for stable.fast algo in pcalg
m.max = 4      #max. number of conditioning variables 
PC_file_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Caulsal_Network_Results/"
length(tissues)
for(i in 12:26){
  tissue_type_to_study = tissues[i]
  print(paste0("The present tissue is:",tissue_type_to_study))
  resid_mat = get_Residualized_GTEx(TPM_thres,low_var_thres,tissue_type_to_study,resid_file_prefix)
  Glasso_obj = get_Glasso_Graph(resid_mat,n.cores,tissue_type_to_study,Glasso_file_prefix)
  get_Bayesian_Network(alpha,numCores,m.max,resid_mat,Glasso_obj,PC_file_prefix)
}
