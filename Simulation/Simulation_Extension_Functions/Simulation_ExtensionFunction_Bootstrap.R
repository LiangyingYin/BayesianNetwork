#########################################################################
# This is used for bootstrap PC algorithm network.
# INPUT:
# samples_data, the dataset we used for learning causal graph
# times_bootstrap, times of bootstrap settings
# cut_off,used for setting if an edge should be remained in the final output.
# 
# OUTPUT: the learned causal network after bootstrap. 
# NOTE: We only remain the edge that can be learned from more than a specific percentage of graph. 
# The percentage is set by 'cut_off' parameters

#########################################################################

bootstrap_pc_network<-function(samples_data,times_bootstrap,cut_off) 
{
  adj_list = list()
  adj_final = NULL
  for(i in 1:times_bootstrap)
  {
    new_sample_index = sample(seq(1,nrow(samples_data)), nrow(samples_data),replace = TRUE)
    new_samples = samples_data[new_sample_index,]
    pc_obj = PC_Algorithm_GeneNetwork(new_samples)
    
    adj = as(pc_obj@graph,"matrix")
    adj_list[[i]] = adj
    if(i==1){adj_final = adj}
    else{ adj_final = adj_final + adj }
    
    pcSimple.fit = GeneToPhenotype_Network(new_samples)
    as.numeric(pcSimple.fit$G)
    
  }
  
  adj_final_output = adj_final
  non_zero_location = which(adj_final!=0,arr.ind = T)
  rows = dim(non_zero_location)[1] #num of row
  columns = dim(non_zero_location)[2] #num of column
  for (i in 1:rows)
  {
    x_loc = non_zero_location[i,1]
    y_loc = non_zero_location[i,2]
    edge_mark = adj_final[x_loc,y_loc]
    if(edge_mark>=times_bootstrap*cut_off)
    {
      adj_final_output[x_loc,y_loc] = 1
    }
    else { adj_final_output[x_loc,y_loc] = 0 }
  }
  
  graph_final_est = as(adj_final_output,"graphNEL")
  return(graph_final_est)
}


