
#***************************************************************************************
# Written by: yinly, August 6,2020
# Description: Calculate total and direct effect of genes on the studied outcome based
# on our estimated causal network that only includes nodes somehow connected to outcome
#***************************************************************************************

library(pcalg)
library(glassoFast)
library(glasso)
library(dplyr)
library(igraph)
library(graph)
library(Matrix)

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
  estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov_mat,graphEst=graph.fit,technique="RRC",type="cpdag")
  #estimate_jointIDA_list = ida(x.pos=x.pos, y.pos=y.pos, cov_mat, graph.fit, method = c("optimal"),
  #                             y.notparent = FALSE, verbose = FALSE, all.dags = NA, type = c("cpdag") )
  # Iterate all candidate genes and calculate their IDA and jointIDA with regard to the studied phenotype
  
  # Create a matrix to store the IDA and jointIDA for each causal genes
  #Effects_mat = matrix(0,nrow=length(x.pos),ncol=3)
  Effects_mat = matrix(0,nrow=length(x.pos),ncol=2)
  
  for(i in 1:length(x.pos)){
    node1 = nodes_estimateGraph[i] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    IDA.estimate <- ida(node1,node2, cov_mat, graph.fit, method = "local", type = "cpdag")
    #IDA.estimate <- idaFast(node1,node2, cov_mat,graph.fit)	
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

Get_Causal_effects_of_outcome_related_genes <- function(bgadd,pc_graph_name,Selected_genes_name,target_outcome_name,Merged_graph_name){
  #*************************************************************************************************
  # Step 1: load merged graph and select nodes that directly or indirectly connect to the outcome
  #*************************************************************************************************
  
  # Get the merged graph for studied trait
  estimateMatrix = bgadd  
  estimateMatrix = t(estimateMatrix) ## to make bgadd consistent with the graph derived from GTEx
  
  # Keep the names of the last column and row consistent, i.e., 71 as the name for "outcome"
  rownames(estimateMatrix) = seq(1,nrow(estimateMatrix))
  colnames(estimateMatrix) = seq(1,nrow(estimateMatrix))
  
  # Convert the combined graph as a igraph object
  graph.fit = as(estimateMatrix, "graphNEL")
  igraph_fit = graph_from_graphnel(graph.fit, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
  #igraph_fit = as.igraph(graph.fit)

  # Selected nodes directly or indirectly connected to the outcome
  # Get the distance table for the igraph.fit
  igraph_fit_distance = distances(igraph_fit)
  outcome_distance = igraph_fit_distance[,ncol(igraph_fit_distance)]
  selected_index = which(outcome_distance!=Inf)
  # Add outcome column to the selected columns for the graph
  selected_gene_index = selected_index[-length(selected_index)]

  # Select nodes that somehow connected to the outcome column
  selectedMatrix = estimateMatrix[selected_index,selected_index]
  rownames(selectedMatrix) = seq(1,length(selected_index))
  colnames(selectedMatrix) = seq(1,length(selected_index))
  cat("There are",nrow(selectedMatrix),"nodes remaining after distance filteration for",target_outcome_name,"\n")
  # Convert the graph with selected nodes as a igraph object
  graph.fit = as(selectedMatrix, "graphNEL")
  
  #convert format of the graph that have negative edges
  estimateGraph = graph.adjacency(selectedMatrix, mode="directed", weighted=TRUE)
  estimateGraph = igraph.to.graphNEL(estimateGraph)
  #**************************************************************************
  # Step 2: calculate covariance matrix for the graph with selected nodes
  #**************************************************************************
  # load Selected_genes with zMin for UKBB graph
  Selected_genes_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Merged_Genes_Results_Update/"
  #Selected_genes_name = "Whole_Blood_CAD_Selected_Genes.csv"
  Selected_genes_all = fread(paste0(Selected_genes_prefix,Selected_genes_name))
  # pick out genes that directly or indirectly connected to the outcome based on distance_table 
  Selected_genes = Selected_genes_all[selected_gene_index,]
  
  #load imputed gene expression data and real phenotypes from UKBB
  expression_file_dir = "/exeh_3/yinly/BayesianNetwork/01_PrediXican/Result/Filtered/UKBB_Whole_Blood_CIS_Trans_expression.txt"
  expr = fread(expression_file_dir )
  # load genetically imputed clinical trait
  predicted_pheno_dir = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/02_Multiple_Traits/Genetically_Predicted_BMI_HbA1c.csv"
  predicted_pheno = fread(predicted_pheno_dir)
  #colnames(predicted_pheno) = c("FID","ENSG_LDL") # The second element value needs to be revised accordingly
  #colnames(predicted_pheno) = c("FID",predicted_trait_name)
  expr_new = inner_join(expr, predicted_pheno,by="FID")
  
  # pick out selected genes
  selected_genes_Ensem = Selected_genes$genes_names
  selected_genes_colnames = c("FID",selected_genes_Ensem)
  expr_selected = subset(expr_new,select = selected_genes_colnames)

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
  # Calculate covariance matrix based on selected genes
  Cov_Estimated = cov(UKBB_Valid)
  
  # Store the calculated covariance matrix for studied tissue-trait pairs
  # Covariance_Prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Covariances/"
  # Covariance_filename = gsub("_Selected_Genes.csv","",Selected_genes_name)
  # fwrite(as.data.frame(Cov_Estimated),paste0(Covariance_Prefix,Covariance_filename,"_Covariance.csv"))
  
  #************************************************************************************
  # Step 3: Adjust covariance matrix based on pMax matrix  
  #************************************************************************************
  pheno_num = 1
  gene_num =  nrow(Selected_genes)
  rho_pMax_len = gene_num + pheno_num
  rho_pMax = matrix(0,nrow= rho_pMax_len,ncol=rho_pMax_len)
  
  # load pc_graph(GTEx causal graph with pMaxs)
  # note: pc_graph is a list, the second component is the pMax matrix
  #pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_Network_Relax_Rand/"
  #pc_graph_name = "Whole_Blood_CAD_Causal_Network_PCSelect.Rdata"
  load(paste0(pc_graph_prefix,pc_graph_name))
  
  #***********************************************************************************
  # Revised by yinly, Oct 12, 2020
  # Description: update the pMax matrix derived from the pc algorithm accordingly
  #***********************************************************************************
  pMax_temp = pc_obj@pMax
  pMax_temp = forceSymmetric(pMax_temp,uplo="U")
  pMax_update = pMax_temp
  pMax_update[which(pMax_update==-1)] = 1
  adj_temp = as(pc_obj@graph,"matrix")
  alpha = 0.05
  pMax_update[which(pMax_update<alpha & adj_temp ==0)] = 1  
  #pMax_final = Matrix(pMax_update,nrow=nrow(pMax_update))
  pMax_final = as.matrix(pMax_update)
  pc_obj@pMax = pMax_final
  # end of revision
  
  # Get rho_pMax matrix for the overall causal graph
  #rho_pMax_gene = pc_graph[[2]] ##The 2nd element of pc_graph is the pMax matrix
  rho_pMax_gene = pc_obj@pMax  ## pc_obj is directly derived from pc algorithm, there is a slot named pMax
  #***********************************************************************************
  # Add pMax row and column for 2 genetically predicted clinical traits
  #***********************************************************************************
  # Add the 1st genetically predicted clinical trait
  zMin_trait1 = Selected_genes_all$zMin.x[1:(nrow(Selected_genes_all)-1)]
  rho_pMax_predicted_trait1 = 2*pnorm(-abs(zMin_trait1))
  rho_pMax_predicted_trait1[which(is.na(rho_pMax_predicted_trait1))] = -1
  len_predictors = length(rho_pMax_predicted_trait1)
  rho_pMax_predictors = matrix(0,nrow= len_predictors,ncol=len_predictors)
  rho_pMax_predictors[1:(len_predictors-1),1:(len_predictors-1)] = rho_pMax_gene
  rho_pMax_predictors[len_predictors,] = rho_pMax_predicted_trait1
  rho_pMax_predictors[,len_predictors] = rho_pMax_predicted_trait1

  # Add the 2nd genetically predicted clinical trait
  zMin_trait2 = Selected_genes_all$zMin.y[1:nrow(Selected_genes_all)]
  rho_pMax_predicted_trait2 = 2*pnorm(-abs(zMin_trait2))
  rho_pMax_predicted_trait2[which(is.na(rho_pMax_predicted_trait2))] = -1
  len_predictors2 = length(rho_pMax_predicted_trait2)
  rho_pMax_predictors2 = matrix(0,nrow= len_predictors2,ncol=len_predictors2)
  rho_pMax_predictors2[1:(len_predictors2-1),1:(len_predictors2-1)] = rho_pMax_predictors
  rho_pMax_predictors2[len_predictors2,] = rho_pMax_predicted_trait2
  rho_pMax_predictors2[,len_predictors2] = rho_pMax_predicted_trait2
  #diag(rho_pMax_predictors2) <- 0

  # pick out the columns that are directly or indirectly connected to the outcome column
  rho_pMax_gene = rho_pMax_predictors2[selected_gene_index,selected_gene_index]
  rho_pMax_gene_outcome = 2*pnorm(-abs(Selected_genes$zMin.z)) ##
  rho_pMax_gene_outcome[which(is.na(rho_pMax_gene_outcome))] = -1
  rho_pMax[1:gene_num,1:gene_num] = rho_pMax_gene
  rho_pMax[,rho_pMax_len] = c(rho_pMax_gene_outcome,-1)
  rho_pMax[rho_pMax_len,] = c(rho_pMax_gene_outcome,-1)
  diag(rho_pMax) <- 0
  rho_pMax[rho_pMax==-1] <- 1

  # Adjust covariance matrix based on the pMax matrix
  Cov_Estimated_adjusted = glassoFast(Cov_Estimated, rho=rho_pMax)$w

  # File prefix for causal effects files
  Effects_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Results_Effects_Update/"
  Effects_filename = gsub(".Rdata","",Merged_graph_name)
  
  # Calculated IDA and jointIDA based on pvalue-adjusted covariance matrix
  Effects_mat_Pvalue = NULL
  Effects_mat_Pvalue = Get_IDA_and_jointIDA(estimateGraph,Cov_Estimated_adjusted,graph.fit,selected_genes_Ensem)
  
  Effects_pvalue_filepath = paste0(Effects_prefix,Effects_filename,"_Pvalue.csv")
  fwrite(as.data.frame(Effects_mat_Pvalue),Effects_pvalue_filepath)
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
  
  # Calculate IDA and jointIDA based on iterative pvalue-adjusted covariance matrix
  Effects_mat_Pvalue_Iterative = NULL
  Effects_mat_Pvalue_Iterative = Get_IDA_and_jointIDA(estimateGraph,Cov_Estimated_pvalue_iterative_adjusted,graph.fit,selected_genes_Ensem)
  
  Effects_iterative_pvalue_filepath = paste0(Effects_prefix,Effects_filename,"_Iterative_Pvalue.csv")
  fwrite(as.data.frame(Effects_mat_Pvalue_Iterative),Effects_iterative_pvalue_filepath)
  
}

#Merged_graph_prefix = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/02_Multiple_Traits/Merged_Graphs/"
#Merged_graph_name = "LDL_CAD_Whole_Blood.Rdata"
#load(paste0(Merged_graph_prefix,Merged_graph_name)) # loaded object name bgadd2

#pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_Network_Relax_Rand/"
#pc_graph_name = "Whole_Blood_LDL_Causal_Network_PCSelect.Rdata"
#Selected_genes_name = "UKBB_LDL_CAD_Whole_Blood.csv"
#target_outcome_name="30780-0.0" # for DM(2443-0.0), LDL(30780-0.0)
#Get_Causal_effects_of_outcome_related_genes(bgadd,pc_graph_name,Selected_genes_name,target_outcome_name,Merged_graph_name)

Merged_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Merged_Graphs_Update/"
pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_Network_Update_Valid/"
Info_filepath = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Mutiple_Traits/Effects_Calculation_Dict.csv"
Info_mat = as.matrix(fread(Info_filepath))
#for(i in 1:nrow(Info_mat)){
for(i in 6:7){
     Merged_graph_name = Info_mat[i,1]
     load(paste0(Merged_graph_prefix,Merged_graph_name)) # loaded object name bgadd2
     pc_graph_name = Info_mat[i,2]
     Selected_genes_name = Info_mat[i,3]
     target_outcome_name = Info_mat[i,4]
     predicted_trait_name = Info_mat[i,5]
     Get_Causal_effects_of_outcome_related_genes(bgadd,pc_graph_name,Selected_genes_name,target_outcome_name,Merged_graph_name,predicted_trait_name)
}
#Merged_graph_name = "BMI_HbA1c_COVID_Confirm_Whole_Blood.Rdata"
#load(paste0(Merged_graph_prefix,Merged_graph_name)) # loaded object name bgadd
#pc_graph_name = "Whole_Blood_AF_Causal_Network_PCSelect.Rdata"
#Selected_genes_name = "UKBB_BMI_HbA1c_COVID_Confirm_Whole_Blood.csv"
#target_outcome_name = "COVID_Confirm"
#predicted_trait_name = c("ENSG_BMI","ENSG_Diabetes")
#Get_Causal_effects_of_outcome_related_genes(bgadd,pc_graph_name,Selected_genes_name,target_outcome_name,Merged_graph_name)
