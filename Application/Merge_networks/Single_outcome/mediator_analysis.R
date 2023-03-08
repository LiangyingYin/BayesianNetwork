library(pcalg)
library(CePa)
library(reshape2)
library(data.table)
library(dplyr)
library(coop)
library(parallel)
library(graph)
library(grid)
library(pcalg)
library(Rgraphviz)
library(Corbi)
library(glmnet)
library(MLmetrics)
library(MASS)
library(igraph)  
library(stats)
library(RNOmni)
library(ParallelPC)
library(glasso)
library(glassoFast)
library(RBGL)
library(epiR)
library(RiskPortfolios)
library(ParallelPC)
library(glasso)
library(glassoFast)
library(RBGL)
library(epiR)
library(matrixcalc)
library(CVglasso)
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(dplyr)
library(coop)
library(pcalg)
library(parallel)
library(glasso)
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(pcalg)
library(coop)
library(parallel)
library(rlist)

library(pcalg)
library(CePa)
library(reshape2)
library(data.table)
library(dplyr)
library(coop)
library(parallel)
library(graph)
library(grid)
library(pcalg)
library(Rgraphviz)
library(Corbi)
library(glmnet)
library(MLmetrics)
source("/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Exploratory_Analysis/Mediator_genes_detection.R")

Merged_graph_prefix <- "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Merged_Graphs_Update/"
Selected_genes_prefix <- "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_Selected_Genes_Update/"
Mediator_analysis_prefix <- "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Mediation_Analysis_Update/"
Info_filepath <- "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/GTEx_UKBB_Merge_Dict_Update_1April.csv"
Info_mat <- as.matrix(fread(Info_filepath))

#index <- c(12,24,36,67,125)
index <- c(53,54,86,87)
#for(i in 1:nrow(Info_mat)){
for(j in 1:length(index)){
  i = index[j]
  #i <- 95
  Merged_graph_name <- Info_mat[i,5] ##the column index for merged graph
  # load bgadd2(whole causal graph)
  #Merged_graph_name = "Whole_Blood_CAD_Graph.Rdata"
  load(paste0(Merged_graph_prefix,Merged_graph_name))##loaded object name: bgadd2
  Selected_genes_name <- Info_mat[i,4] 
  #*************************************************************************************************
  # Step 1: load merged graph and select nodes that directly or indirectly connect to the outcome
  #*************************************************************************************************
  
  # Get the merged graph for studied trait
  estimateMatrix <- bgadd2  
  estimateMatrix <- t(estimateMatrix) ## to make bgadd2 consistent with the graph derived from GTEx
  
  # Keep the names of the last column and row consistent, i.e., 71 as the name for "outcome"
  rownames(estimateMatrix) <- seq(1,nrow(estimateMatrix))
  colnames(estimateMatrix) <- seq(1,nrow(estimateMatrix))
  adj <- estimateMatrix
  mediator_result <- Check_Mediator2(adj) 

  # load Selected_genes with zMin for UKBB graph
  #Selected_genes_name = "Whole_Blood_CAD_Selected_Genes.csv"
  Selected_genes_all <- fread(paste0(Selected_genes_prefix,Selected_genes_name))
  # save mediator detection results to a matrix
  ifmediator_gene <- as.vector(unlist(mediator_result[[1]]))
  gene_type <- as.vector(unlist(mediator_result[[2]]))
  detection_result <- cbind(Selected_genes_all,ifmediator_gene,gene_type)
  mediator_filename = gsub("Selected_Genes","Mediator_Analysis",Selected_genes_name)
  fwrite(as.data.frame(detection_result),paste0(Mediator_analysis_prefix,mediator_filename))
}
