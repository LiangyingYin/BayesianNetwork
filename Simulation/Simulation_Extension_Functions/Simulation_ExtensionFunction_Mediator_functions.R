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

###########################Begin######################################################
# VersionV1.0 added by fyn 29/3/2021
#Make the adj to 1 or 0, without the weight
Delete_Weight_For_adj<-function(adj)
{
  #temp comment: edge with weight to be 1
  adj2 = adj
  non_zero_location2 = which(adj!=0,arr.ind = T)
  for(i in 1:dim(non_zero_location2)[1])
  {
    a = non_zero_location2[i,1] # no-zero location in row
    b = non_zero_location2[i,2] # no-zero location in column
    adj2[a,b] = 1
  }
  return(adj2)
}
###########################End########################################################

###########################Begin######################################################
###Check all the node in a graph, if each of them has mediator refer to outcome#######
Check_Mediator <- function(adj)
{
  adj_gene_outcome = Delete_Weight_For_adj(adj)
  graph.fit = as(adj_gene_outcome, "graphNEL")
  igraph_fit = graph_from_graphnel(graph.fit, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
  igraph_fit_distance_true = distances(igraph_fit)
  distance_list_true = igraph_fit_distance_true[,ncol(igraph_fit_distance_true)]
  
  # to do: need to know the meaning of the output results
  igraph_fit_distance = distances(igraph_fit)
  outcome_distance = igraph_fit_distance[1:(nrow(igraph_fit_distance)-1),ncol(igraph_fit_distance)]
  selected_index = which(outcome_distance!=Inf)
  
  true_distance_list = outcome_distance
  true_adj_list = adj[1:(nrow(adj)-1),ncol(adj)]
  
  for(i in 1:length(true_adj_list))
  {
    weight = true_adj_list[i]
    if(weight==0)
    {
      true_adj_list[i]=0
    }
    else
    {
      true_adj_list[i]=1
    }
  }
  
  differ_list = true_distance_list-true_adj_list
  ifmediator_gene = differ_list # for the differ=0, we need to further analyze if it's mediator or not.the value is the mediator num
  gene_type = differ_list # gene_type, 2 for both direct&indirect, 0 for only mediator, 1 for only direct. to outcome, 3 for non
  
  for(i in 1:length(differ_list))
  {
    ifmediator = differ_list[i]
    #if its zero,we need to make sure if it only have direct causal relation, or both for direct and undirect causal relation.
    if(ifmediator==0)
    {
      adj_tmp = adj
      adj_tmp[i,ncol(adj_tmp)] = 0
     
      
      igraph = graph.adjacency(adj_tmp, mode="directed")
      path = get.shortest.paths(igraph,i,ncol(adj_tmp))$vpath[[1]]
      if(length(path)==0) # then there is no other path between this gene to the outcome.
      {
        ifmediator_gene[i]=0
        gene_type[i]=1 # only direct relation
      }
      else# there is connection between this gene to outcome by mediator
      {
        ifmediator_gene[i] = length(path)-1
        gene_type[i]=2
      }
    }
    else if(ifmediator==Inf) # for these part, we need to make sure if its Inf or not
    {
      ifmediator_gene[i]=0
      gene_type[i]=3 # for mediator, this gene is related to outcome only by mediator
    }
    else
    {
      gene_type[i]=0 # these genes are related to outcome only by mediator.
    }
  }
  
  result_list = list()
  result_list[[1]] = ifmediator_gene
  result_list[[2]] = gene_type
  return(result_list)
}

Check_Mediator2 <- function(adj)
{
  adj_gene_outcome = Delete_Weight_For_adj(adj)
  # adj_tmp = adj_gene_outcome
  ifmediator_gene = seq(1,nrow(adj)) # for the differ=0, we need to further analyze if it's mediator or not.the value is the mediator num
  gene_type = seq(1,nrow(adj)) # gene_type, 2 for both direct&indirect, 0 for only mediator, 1 for only direct. to outcome, 3 for non
  ifmediator_gene[]=0
  gene_type[]=0
  
  for(i in 1:nrow(adj_gene_outcome))
  {
    adj_tmp = adj_gene_outcome
    if(adj_tmp[i,ncol(adj_gene_outcome)]==0)
    {
      igraph = graph.adjacency(adj_tmp, mode="directed")
      path = get.shortest.paths(igraph,i,ncol(adj_tmp))$vpath[[1]]
      if(length(path)==0)# which shows that there is no path between this gene to outcome
      {
        ifmediator_gene[i] = 0
        gene_type[i]=3 # no any relation to outcome
      }
      else if(length(path)!=0)
      {
        ifmediator_gene[i] = 1 # if you want to count the mediator, then maybe should change the function we use
        gene_type[i]= 0 # only have mediator to outcome
      }
    }  
    else if(adj_tmp[i,ncol(adj_gene_outcome)]==1)
    {
      adj_tmp[i,ncol(adj_gene_outcome)] = 0 # delete the direct casual relation.
      igraph = graph.adjacency(adj_tmp, mode="directed")
      path = get.shortest.paths(igraph,i,ncol(adj_tmp))$vpath[[1]]
      if(length(path)==0)# which shows that there is no path between this gene to outcome
      {
        ifmediator_gene[i] = 0
        gene_type[i]=1 # only direct
      }
      else if(length(path)!=0)
      {
        ifmediator_gene[i] = 1
        gene_type[i]=2  #both direct&indirect
      }
    }
    
  }
  result_list = list()
  result_list[[1]] = ifmediator_gene
  result_list[[2]] = gene_type
  return(result_list)
}
###########################End########################################################

###########################Begin######################################################
#Standardize the mediator results, ingore their mediator type, Ture-have mediator, FALSE-dont have mediator.
Standardize_Mediator_result<-function(ifmediator_gene)
{
  ifmediator_gene_std = ifmediator_gene
  for(i in 1:length(ifmediator_gene))
  {
    if(ifmediator_gene[i]==0)
    {
      ifmediator_gene_std[i] = FALSE
    }
    else
    {
      ifmediator_gene_std[i] = TRUE
    }
  }
  return(ifmediator_gene_std)
}
###########################End########################################################