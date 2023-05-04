

#**************************************************************
# Organized by yinly, May 12,2020
#**************************************************************

#**********************************************************************************************
# Goal: Merge GTEx gene network and UKBB gene-phenotype network
# pc_graph: object derived from bootstrapped pc([[1]]adjacency matrix,[[2]]pMax matrix)
# pcSimple.fit: object derived from pcSimple_Select
#**********************************************************************************************
MergeGetx_UKBB<-function(pc_graph,pcSimple.fit,Selected_genes)
{
  #*************************************
  # This function adds background knowledge to a graph from pcfit
  #****************************************
  library(graph)
  library(pcalg)
  
  #*******************************************************************
  # get adjacency matrix for GTEx network derived from bootstrapped pc
  #*******************************************************************
  # pc_graph = pc_graph@graph 
  pc_graph_adj = pc_graph[[1]]
  adj_mat = as(pc_graph_adj, "matrix") ##[i,j] =1 means an edge from i to j (graph package documentation)
  
  #*******************************************************************************
  # get adjacency matrix for UKBB(Select top N genes based on simulation results)
  #*******************************************************************************
  # Select the top K if the pcSimple.fit has more than K edgeds
  # Get selected genes from pcSimple.fit
  # genes_names = names(pcSimple.fit$G)
  # assoc =  as.vector(unlist(pcSimple.fit$G))
  # zMin = as.vector(unlist(pcSimple.fit$zMin))
  # tissue_trait_result = data.frame(genes_names,assoc,zMin,stringsAsFactors=FALSE)
  # Genes_result_selected = tissue_trait_result[tissue_trait_result$assoc==TRUE,]
  # Genes_selected_sort = Genes_result_selected[order(-Genes_result_selected$zMin),]
  # if(nrow(Genes_selected_sort)>70){
  #   gene_to_outcome = Genes_selected_sort[1:70,2]
  # }else{
  #   gene_to_outcome = Genes_selected_sort[,2]
  # }
  gene_to_outcome = as.numeric(Selected_genes$assoc)
  #gene_to_outcome = as.numeric(pcSimple.fit$G) ##Further revision is required based on simulation
  
  #*******************************************************************************************
  # Commented by:yinly: when combine adjacency matrix derived from GTEx with that from UKBB, 
  # the gene_to_outcome was set as the first column, this may affect the finally derived graph
  #*******************************************************************************************
  # final_mat = cbind(gene_to_outcome,adj_mat) ##[i,j] means an edge from i to j (here it is gene -> disease, but NOT the other way round
  # final_mat <- rbind(rep(0,ncol(adj_mat)+1),final_mat) ##need to check if the diagonals are all zero 
  
  # Revised by yinly: adjacency matrix derived from GTEx was combined first
  final_mat = cbind(adj_mat,gene_to_outcome) ##[i,j] means an edge from i to j (here it is gene -> disease, but NOT the other way round
  final_mat <- rbind(final_mat,rep(0,ncol(adj_mat)+1)) ##need to check if the diagonals are all zero 
  #pc_graph_adj2 = as(adj_mat2, "graphNEL")
  
  
  #*******************************************************
  # Merge two graph by add background knowledge
  #*******************************************************
  #the gInput takes [i,j] to mean an edge from j to i, which is opposite to the convention in graph package
  mat_input = t(final_mat) 
  
  t1 = proc.time()
  bgadd2 = addBgKnowledge(gInput = mat_input,  
                          verbose = TRUE, 
                          checkInput = FALSE)					   
  
  proc.time()-t1
  
  return(bgadd2)
}

library(data.table)
GTEx_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Reoriented_GTEx_Network_Update/"
UKBB_prefix = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/Order3_Results_Update/"
Selected_genes_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_Selected_Genes/"
# The dictionary that match the corresponding GTEx gene-gene network with UKBB gene-outcome
file_Dict = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/GTEx_UKBB_Merge_Dict.csv"
Dict = as.matrix(fread(file_Dict))
for(i in 1:nrow(Dict)){
   # load GTEx object derived bootstrapped pc
   #i = 105 
   GTEx_name = Dict[i,2]
   #GTEx_name = "Whole_Blood_CAD_Causal_Network_PCSelect.Rdata"
   load(paste0(GTEx_prefix,GTEx_name))##object name: pc_graph
   # load pcSimple.fit
   UKBB_name = Dict[i,3]
   #UKBB_name = "UKBB_CAD.Rdata"
   load(paste0(UKBB_prefix,UKBB_name))##object name: pcSimple.fit
   Selected_genes_name = Dict[i,4]
   #Selected_genes_name = "Whole_Blood_CAD_Selected_Genes.csv"
   Selected_genes = fread(paste0(Selected_genes_prefix,Selected_genes_name))
   # Merge GTEx object with pcSimple.fit based on selected genes
   bgadd2 = MergeGetx_UKBB(pc_graph,pcSimple.fit,Selected_genes)
   Merged_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Merged_Graphs/"
   graph_name = paste0(gsub("Causal_Network_PCSelect.Rdata","",GTEx_name),"Graph.Rdata")
   save(bgadd2,file=paste0(Merged_graph_prefix,graph_name))
}

#*****************************************************************************
# The following code is for testing only
#*****************************************************************************
#GTEx_name="Artery_Coronary_CAD_Causal_Network_PCSelect.Rdata"
#load(paste0(GTEx_prefix,GTEx_name))##object name: pc_graph
#UKBB_name = "UKBB_CAD_Artery_Coronary.Rdata"
#load(paste0(UKBB_prefix,UKBB_name))##object name: pcSimple.fit
#Selected_genes_name = "Artery_Coronary_CAD_Selected_Genes.csv"
#Selected_genes = fread(paste0(Selected_genes_prefix,Selected_genes_name))
# Merge GTEx object with pcSimple.fit based on selected genes
#Merged_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Results_Graphs_New/"
#bgadd2 = MergeGetx_UKBB(pc_graph,pcSimple.fit,Selected_genes)
#graph_name = paste0(gsub("Causal_Network_PCSelect.Rdata","",GTEx_name),"Graph.Rdata")
#save(bgadd2,file=paste0(Merged_graph_prefix,graph_name))
