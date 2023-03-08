
#**********************************************************************************************************
#  Written by: yinly, June 7, 2020
#  Description: Exploratory analysis based on PCSelect
#  Revised by yinly, June 9, 2021
#  Description: for each tissue, we only need to infer one gene-gene network!
#**********************************************************************************************************


#***************************************
#  pc algorithm applied to GTEx 
#***************************************
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(dplyr)
library(coop) 
library(pcalg)
library(parallel) 
library(glasso)
library(qgraph)


setwd("/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/Residmat_All_Update/")
pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_Network_Update/"
files = list.files(path = ".")
#index = c(1,12,23,24,25,34,36,39,50,52,54,65,66)
#index = c(13,33)
index <- c(55,56,89,90)
#for(i in 2:4){
for(j in 1:length(index)){
    i <- index[j]
    #i <- 1
    # Load input data resid_mat_selected
    load(files[i])
    # Apply glassoFast on resid_mat_selected to remove unlikely edges
    n.cores  = 15 
    cov_mat = coop::covar(resid_mat_selected)
    rho_theoretical =  sqrt( log(ncol(resid_mat_selected)) / nrow(resid_mat_selected) ) 
    print(rho_theoretical)		
    library(glassoFast)	
    glasso_obj = glassoFast(S=cov_mat, rho= rho_theoretical) 
    # Use PC algorithm to infer GTEx network
    corr_mat = pcor(resid_mat_selected) 
    suffStat <- list(C=corr_mat ,n=nrow(resid_mat_selected)) ##pcor is a function from "coop" for fast computation of correlation matrix
    #**********************
    # enforce symmetricity in the partial correlation matrix obtained from glasso
    #**********************************
    library(Matrix)
    FixedGaps_mat = forceSymmetric(glasso_obj$wi)
    #https://stackoverflow.com/questions/18165320/creating-a-symmetric-matrix-in-r

    FixedGaps_mat[FixedGaps_mat!=0] <- 999 
    FixedGaps_mat[FixedGaps_mat==0] <- 1
    FixedGaps_mat[FixedGaps_mat==999] <- 0
    FixedGaps_mat = matrix( as.logical(FixedGaps_mat), nrow = nrow(FixedGaps_mat) ) 

    t1=proc.time()
    alpha = 0.05   #alpha (p-value) threshold for 
    numCores = 20  #number of cores for stable.fast algo in pcalg
    m.max = 4      #max. number of conditioning variables 
    pc_obj = pc(suffStat, 
                indepTest=gaussCItest, 
                p=ncol(corr_mat), 
                skel.method="stable.fast",
                fixedGaps = FixedGaps_mat,
                u2pd = "relaxed",
                maj.rule = TRUE,
                solve.confl = TRUE,
                alpha = alpha, 
                numCores = numCores ,
                m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209  
                verbose= FALSE)
    proc.time()-t1
    file_name = gsub("_resid_mat_selected.Rdata","",files[i])
    save(pc_obj, file = paste0(pc_graph_prefix,file_name,"_Causal_Network_PCSelect.Rdata")) 
    cat("The inference of the",i,"PC graph",file_name,"is completed!","\n")   
    file_name = gsub("_Causal_Network_PCSelect.Rdata","",files[i])
    adj_temp = as(pc_obj@graph,"matrix")
    amat = t(adj_temp)
    cat("This is the",i,"graph",file_name,"\n")
    isValidGraph(amat,verbose=TRUE)
}


