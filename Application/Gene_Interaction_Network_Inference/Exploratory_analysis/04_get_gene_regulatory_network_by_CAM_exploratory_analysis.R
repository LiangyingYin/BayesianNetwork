
library(data.table)
library(CAM)

setwd("/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/")
Resid_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/Residmat_All_Update/"
CAM_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/CAM_Network_Update/" # directory stored the inferred gene interaction network by causal additive model
invalid_Info = as.matrix(fread("InvalidGraphs_Update.csv.csv"))
for(i in 1:nrow(invalid_Info)){
    #invalid_resid_name = invalid_Info[i,2]
    #i <- 3
    #i <- 3
    load(paste0(Resid_prefix,invalid_Info[i,5]))
    cam_obj = CAM(resid_mat_selected,scoreName = "SEMLIN",parsScore = list(numBasisFcts=1), variableSel = TRUE, variableSelMethod = selLasso,pruning= TRUE, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts=1),numCores = 20)
    file_name = gsub("_resid_mat_selected.Rdata","",invalid_Info[i,5])
    save(cam_obj, file = paste0(CAM_prefix,file_name,"_CAM.Rdata")) 
    cat("The inference of the",i,"CAM graph",file_name,"is completed!","\n")
}
