#############################################################################
# Merge Gtex and UKBB data
#############################################################################
MergeGetx_UKBB<-function(pcSimple.fit,graph_est)
{
  #*************************************
  # This function adds background knowledge to a graph from pcfit
  #****************************************
  library(graph)
  library(pcalg)
  
 
  outcome = "BMI"   ##name of outcome to be specified here
  
  pc_graph = graph_est
  adj_mat = as(pc_graph, "matrix")  ##[i,j] =1 means an edge from i to j (graph package documentation)
  
  #**********************************
  #   load gene -> outcome links 
  #**********************************		
  # load("pcSimple_Stable_fit_UKBB.Rdata")
  
  #****************************************
  # to do: add to the orig. adjacency matrix (of expressions) gene->disease connections
  
  # outcome G1  G2  G3 ...
  # outcome
  # G1
  # G2
  # G3
  
  gene_to_outcome = as.numeric(pcSimple.fit$G) 
  final_mat = cbind(gene_to_outcome, adj_mat)   ##[i,j] means an edge from i to j (here it is gene -> disease, but NOT the other way round
  final_mat <- rbind(  rep(0,ncol(adj_mat)+1)  , final_mat) ##need to check if the diagonals are all zero 
  
  
  #***********************
  # add background knowledge
  #*********************
  mat_input = t(final_mat)  ##reason for transposing is that the gInput takes [i,j] to mean an edge from j to i, which is opposite to the convention in graph package
  
  t1 = proc.time()
  bgadd2 = addBgKnowledge(gInput = mat_input,  ##reason for transposing is that the gInput takes [i,j] to mean an edge from j to i, which is opposite to the convention in graph package       
                          verbose = TRUE, 
                          checkInput = FALSE)					   
  
  proc.time()-t1
  	
  
  return(bgadd2)
}

