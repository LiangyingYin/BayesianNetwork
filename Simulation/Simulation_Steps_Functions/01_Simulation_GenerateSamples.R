##############################Begin###############################################
# Previous Version for DAG. updated by fyn. 14/5/2021
Simulate_Network<-function(s,g,p,node_num,all_serial,gene_index,snp_index) #prop_direct_causal_genes=prop_direct_causal_genes)
{
  
  rDAG <- randomDAG(n= node_num, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)
  adj = as(rDAG, "matrix") #adj[i,j]!=0 represent that there is an edge i->j

  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  rDAG = igraph.to.graphNEL(rDAG)
  
  tmp_loop = g* prop_direct_causal_genes 

  for (lp in 1:tmp_loop)
  {
    tmp_x = sample(gene_index, size = 1,replace = FALSE)

    tmp_y = n3 #only one phenotype
    
    adj_tmp = as(rDAG, "matrix") 
    
    weight_tmp = runif(1, min = min_weight_gene_to_pheno, max = max_weight_gene_to_pheno)
    adj[tmp_x,tmp_y] = weight_tmp
    
  }
  
  rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
  
  rDAG = igraph.to.graphNEL(rDAG)
  
  non_zero_location = which(adj!=0,arr.ind = T)
  rows = dim(non_zero_location)[1] #num of row
  columns = dim(non_zero_location)[2] #num of column
  
  #pre-processing in data
  #set SNP-SNP relation to 0
  #set gene-SNPs relation to 0
  
  for (i in 1:rows)
  {
    #genes and snps are randomly selected from data instead of choosing them in order as before
    
    
    #set SNP-SNP relation to 0
    if((non_zero_location[i,1] %in% snp_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
    #set gene-SNPs relation to 0
    else if((non_zero_location[i,1] %in% gene_index)&&(non_zero_location[i,2] %in% snp_index))
    {
      adj[non_zero_location[i,1],non_zero_location[i,2]]=0
    }
   
  }
  
  #set pheno-other relation to 0
  adj[n3,]=0
  
  adj[snp_index,n3] = 0  #set snps to phenotype relation to 0
  
  
  #ensure that there is at least one snp related to gene.If there is a gene which it does not have related snps, the imputed gene expression data will be the same in all samples.
  for(gene in gene_index)
  {
    if(length(which(adj[snp_index,gene]!=0))==0)
    {
      snp_loc = sample(snp_index, 1, replace = FALSE)
      weight_tmp = runif(1, min = randomDAG_lB, max = randomDAG_uB)
      adj[snp_loc,gene] = weight_tmp
    }
  }
  
  
  #set snp-gene weaker to be consistent with real situation
  times_to_snpGeneWeight = randomDAG_uB/max_weight_snp_to_gene
  adj[snp_index,gene_index] = adj[snp_index,gene_index]/times_to_snpGeneWeight
  
  # adjust adj to amatType format, to make sure graph is DAG after update
  non_zero_location2 = which(adj!=0,arr.ind = T)
  
  #temp comment: edge with weight to be 1
  adj2 = adj
  for(i in 1:dim(non_zero_location2)[1])
  {
    a = non_zero_location2[i,1] # no-zero location in row
    b = non_zero_location2[i,2] # no-zero location in column
    adj2[a,b] = 1
  }
  
  amat = t(adj2) 
  # in amatType, a[i,j] = 0; a[j,i]=1, represent an edge i->j
  # reference https://www.rdocumentation.org/packages/pcalg/versions/2.6-2/topics/amatType
  
  isDAGresult = isValidGraph(amat = amat, type = "dag",verbose = TRUE)
  isCPDAGresult = isValidGraph(amat = amat, type = "cpdag") ## is a valid CPDAG   completed partially directed acyclic graph 
  isPDAGresult = isValidGraph(amat = amat, type = "pdag") ## is a valid PDAG
  
  if(isDAGresult||isCPDAGresult)
  {
    return(adj) #if the graph is DAG or CPDAG, then return.
  }
}


Simulate_Noise<-function(sample_num2,s,g)
{
  Sigma <- diag(s) #snps error : construct a unit diagonal matrix
  eMat_s <- mvrnorm(sample_num2, mu = rep(0, s), Sigma = Sigma) #produce snp error
  
  eMat_g = NULL
  tmp_num = g
  
  for (g_num in 1:tmp_num)
  {
    eMat_g_tmp = rnorm(sample_num2, mean = 0, sd = eMat_g_tmp_rnorm_sd) # set variance of gene data error manually
    eMat_g =cbind(eMat_g,eMat_g_tmp)
  }
  
  eMat_p <- rnorm(sample_num2, mean = 0, sd = 1) # produce phenotype error
  eMat = cbind(eMat_s,eMat_g,eMat_p)
  
  eMat[,gene_index] = eMat_g
  eMat[,node_num] = eMat_p
  eMat[,snp_index]=eMat_s
  
  return(eMat)
}

