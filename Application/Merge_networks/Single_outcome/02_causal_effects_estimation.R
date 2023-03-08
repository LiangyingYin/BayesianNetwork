
library(data.table)

#***************************************************************************************
# Revised by: yinly, June 19,2020
# Description: Calculate IDA and jointIDA based estimated graph and adjusted covariance
# matrix, the covariance matrix could be adjusted by pvalue-only, pvalue-iterative as well
# as other strategies.
#***************************************************************************************
Get_IDA_and_jointIDA <- function(estimateGraph,cov_mat,graph.fit,selected_genes_Ensem){
  # Calculation of IDA
  nodes_estimateGraph = nodes(estimateGraph)
  phenotype_loc = length(nodes_estimateGraph)  # phenotype is the last node in the Graph
  x.pos = nodes_estimateGraph[!nodes_estimateGraph %in% phenotype_loc]
  y.pos = phenotype_loc
  x.pos = as.numeric(x.pos) # change the data type from character to numeric
  y.pos = as.numeric(y.pos)
  estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov_mat,graphEst=estimateGraph,technique="RRC")
  # Iterate all candidate genes and calculate their IDA and jointIDA with regard to the studied phenotype
  
  # Create a matrix to store the IDA and jointIDA for each causal genes
  #Effects_mat = matrix(0,nrow=length(x.pos),ncol=3)
  Effects_mat = matrix(0,nrow=length(x.pos),ncol=2)
  
  for(i in 1: length(x.pos)){
    node1 = nodes_estimateGraph[i] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    IDA.estimate <- ida(node1,node2, cov_mat, graph.fit, method = "local", type = "pdag")	
    Effects_mat[i,1] = selected_genes_Ensem[i]
    # Revised by: yinly, June 23, 2020
    # Description: the length of the derived IDA.estimate may be larger than 1
    if(length(IDA.estimate) == 1){
      Effects_mat[i,2] = IDA.estimate
    } else{
      IDA.estimate_combine = NULL
      for(j in 1: length(IDA.estimate)){
        if(is.null(IDA.estimate_combine)){
          IDA.estimate_combine = IDA.estimate[j]
        }else{
          IDA.estimate_combine = paste(IDA.estimate_combine,IDA.estimate[j],sep = ",")
        }
      }
      Effects_mat[i,2] = IDA.estimate_combine
    }
    
  }
  #Effects_mat[,3] = estimate_jointIDA_list
  #return(Effects_mat)
  Effects_mat_final = cbind(Effects_mat,estimate_jointIDA_list)
  return(Effects_mat_final)
}



#***************************************************************************************
# Written by: yinly, May 18,2020
# Description: Calculate total and direct effect of genes on the studied outcome based
# on our estimated causal network
#***************************************************************************************
Get_Causal_effects_of_genes <- function(Merged_graph_name,pc_graph_name,Selected_genes_name,target_outcome_name){
  
  library(pcalg)
  library(glassoFast)
  library(glasso)
  library(dplyr)
  library(igraph)
  library(graph)
  
  #*******************************************************
  # Step 1: load graph related variables
  #*******************************************************
  # load bgadd2(whole causal graph)
  Merged_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Merged_Graphs_Update/"
  #Merged_graph_name = "Whole_Blood_CAD_Graph.Rdata"
  load(paste0(Merged_graph_prefix,Merged_graph_name))##loaded object name: bgadd2
  
  # load pc_graph(GTEx causal graph with pMaxs)
  # note: pc_graph is a list, the second component is the pMax matrix
  pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_Selected_Update_Second/"
  #pc_graph_name = "Whole_Blood_CAD_Causal_Network_PCSelect.Rdata"
  load(paste0(pc_graph_prefix,pc_graph_name))
  
  # load Selected_genes with zMin for UKBB graph
  Selected_genes_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_Selected_Genes_Update/"
  #Selected_genes_name = "Whole_Blood_CAD_Selected_Genes.csv"
  Selected_genes = fread(paste0(Selected_genes_prefix,Selected_genes_name))
  
  #**************************************************************************
  # Step 2: load imputed gene expression data and real phenotypes from UKBB
  #**************************************************************************
  #expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Whole_Blood_CIS_Trans_expression.txt"
  expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Lung_predicted_expression.txt"
  expr = fread(expression_file_dir )
  # pick out selected genes
  selected_genes_Ensem = Selected_genes$genes_names
  selected_genes_colnames = c("FID",selected_genes_Ensem)
  expr_selected = subset(expr,select = selected_genes_colnames)
  
  pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/01_Single_Trait/Data/UKBB_Extracted_Traits_31March2022.txt"
  # target name for diabetes
  # target_outcome_name = "CAD"
  pheno = fread(pheno_dir)
  outcome = cbind(pheno[,1], pheno[,..target_outcome_name])
  outcome  = data.frame(na.omit(outcome))
  # test here if FID is the first column
  UKBB_genes_trait = inner_join(expr_selected,outcome, by="FID")
  # Remove the FID column
  UKBB_Valid = UKBB_genes_trait[,-1]
  cat("The sample size for",gsub("_Selected_Genes.csv","",Selected_genes_name),"is",nrow(UKBB_Valid),"\n")
  
  #************************************************************************************
  # Step 3: Get the final pMax matrix for covariance matrix adjustmentand overall graph  
  #************************************************************************************
  # Calculate covariance matrix based on selected genes
  Cov_Estimated = cov(UKBB_Valid)
  
  # Store the calculated covariance matrix for studied tissue-trait pairs
  Covariance_Prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Covariances/"
  Covariance_filename = gsub("_Selected_Genes.csv","",Selected_genes_name)
  #fwrite(as.data.frame(Cov_Estimated),paste0(Covariance_Prefix,Covariance_filename,"_Covariance.csv"))

  pheno_num = 1
  gene_num =  nrow(Selected_genes)
  rho_pMax_len = gene_num + pheno_num
  rho_pMax = matrix(0,nrow= rho_pMax_len,ncol=rho_pMax_len)
  
  # Get rho_pMax matrix for the overall causal graph
  rho_pMax_gene = pc_graph[[2]] ##The 2nd element of pc_graph is the pMax matrix
  # Revised by yinly, June 25,2021
  # Description: re-adjust pmax function based on the adjacency matrix, adj==0 & p<0.01 to p =1
  alpha <- 0.01
  graph_estimate <- pc_graph[[1]] ## The 1st elecment of pc_graph is the inferred graph
  adj_temp <- as(graph_estimate,"matrix")
  rho_pMax_gene[which(rho_pMax_gene<alpha & adj_temp ==0)] = 1  
  # End of revision

  rho_pMax_gene_outcome = 2*pnorm(-abs(Selected_genes$zMin)) ##
  rho_pMax[1:gene_num,1:gene_num] = rho_pMax_gene
  rho_pMax[,rho_pMax_len] = c(rho_pMax_gene_outcome,-1)
  rho_pMax[rho_pMax_len,] = c(rho_pMax_gene_outcome,-1)
  
  diag(rho_pMax) <- 0
  rho_pMax[rho_pMax==-1] <- 1
  Cov_Estimated_adjusted = glassoFast(Cov_Estimated, rho=rho_pMax)$w
  
  #*******************************************************************************************
  # Revised by: yinly, June 19, 2020
  # Description: apply iterative pvalue adjustment strategy on the estimated covariance matrix
  #*******************************************************************************************
  iterateNum = 100
  iterate_error2 = array(0, iterateNum)
  thresh = 0.0001
  Cov_Estimated_pvalue_iterative_adjusted = Cov_Estimated
  for(i in 1:iterateNum)
  {
    GLASSO = glassoFast(Cov_Estimated_pvalue_iterative_adjusted, rho=rho_pMax)
    Cov_Estimated_pvalue_iterative_adjusted = GLASSO$w
    iterate_error2[i] = iterate_error2[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0) *
      log(nrow(UKBB_Valid))/2
    if((i>1)&&(abs(iterate_error2[i]-iterate_error2[i-1])< thresh))
    {
      break
    }
  }
  
  # Get the merged graph for studied trait
  estimateMatrix = bgadd2  
  estimateMatrix = t(estimateMatrix) ## to make bgadd2 consistent with the graph derived from GTEx
  
  # Keep the names of the last column and row consistent, i.e., 71 as the name for "outcome"
  rownames(estimateMatrix) = seq(1,rho_pMax_len)
  colnames(estimateMatrix) = seq(1,rho_pMax_len)
  
  # Convert the combined graph as a grahNEL object
  graph.fit = as(estimateMatrix, "graphNEL")
  
  # If we load bgadd2 from Results_Graphs, then the following 3 lines of codes are required
  # tmp_num = gene_num+1 # obtain some gene and phenotype data from estimateMatrix
  # estimateMatrix = cbind(estimateMatrix[,2:tmp_num],estimateMatrix[,1:1])
  # estimateMatrix = rbind2(estimateMatrix[2:tmp_num,],estimateMatrix[1:1,])
  
  #convert format of the graph that have negative edges
  estimateGraph = graph.adjacency(estimateMatrix, mode="directed", weighted=TRUE)
  estimateGraph = igraph.to.graphNEL(estimateGraph)
  
  #*******************************************************************************
  # Step 5: calcuate the IDA and jointIDA of directly causal genes
  #*******************************************************************************
  # File prefix for causal effects files
  Effects_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Results_Effects_Update/"
  Effects_filename = gsub("_Causal_Network_PCSelect.Rdata","",pc_graph_name)

  # Calculated IDA and jointIDA based on covariance matrix
  Effects_mat = NULL
  Effects_mat = Get_IDA_and_jointIDA(estimateGraph,Cov_Estimated,graph.fit,selected_genes_Ensem)
  #colnames(Effects_mat) = c("EnsembleID","IDA","jointIDA")

  Effects_filepath = paste0(Effects_prefix,Effects_filename,".csv")
  fwrite(as.data.frame(Effects_mat),Effects_filepath)
  
  # Calculated IDA and jointIDA based on pvalue-adjusted covariance matrix
  Effects_mat_Pvalue = NULL
  Effects_mat_Pvalue = Get_IDA_and_jointIDA(estimateGraph,Cov_Estimated_adjusted,graph.fit,selected_genes_Ensem)
  
  Effects_pvalue_filepath = paste0(Effects_prefix,Effects_filename,"_Pvalue.csv")
  fwrite(as.data.frame(Effects_mat_Pvalue),Effects_pvalue_filepath)
  
  # Calculate IDA and jointIDA based on iterative pvalue-adjusted covariance matrix
  Effects_mat_Pvalue_Iterative = NULL
  Effects_mat_Pvalue_Iterative = Get_IDA_and_jointIDA(estimateGraph,Cov_Estimated_pvalue_iterative_adjusted,graph.fit,selected_genes_Ensem)
  
  Effects_iterative_pvalue_filepath = paste0(Effects_prefix,Effects_filename,"_Iterative_Pvalue.csv")
  fwrite(as.data.frame(Effects_mat_Pvalue_Iterative),Effects_iterative_pvalue_filepath)
}

# Load csv file that contains the file names required for the calculation of causal effects
Info_filepath = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/GTEx_UKBB_Merge_Dict_Update_1April.csv"
Info_mat = as.matrix(fread(Info_filepath))
#index = c(86,90,96,98,101,106,107,112,114,117,126)
index <- c(53,54) # 53-54 lung, 86-87 whole blood
#for(i in 135:138){
for(j in 1:length(index)){
  #i = 68
  i = index[j]
  Merged_graph_name = Info_mat[i,5] ##the column index for merged graph
  pc_graph_name = Info_mat[i,2] ##the column index for GTEx graph
  Selected_genes_name = Info_mat[i,4] ##the column index for selected genes
  target_outcome_name = Info_mat[i,6] ##the column index for studied trait in UKBB
  Get_Causal_effects_of_genes(Merged_graph_name,pc_graph_name,Selected_genes_name,target_outcome_name)
}

# The following code are for testing only
# Merged_graph_name = "Whole_Blood_CAD_Graph.Rdata"
# pc_graph_name = "Whole_Blood_CAD_Causal_Network_PCSelect.Rdata"
# Selected_genes_name ="Whole_Blood_CAD_Selected_Genes.csv"
# Get_Causal_effects_of_genes(Merged_graph_name,pc_graph_name,Selected_genes_name)

