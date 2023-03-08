
# Written by yinly, May 11

#
# #***************************************
# #  pc algorithm applied to GTEx
# #***************************************
library(CePa)
library(ParallelPC)
library(reshape2)
library(data.table)
library(dplyr)
library(coop)
library(pcalg)
library(parallel)
library(glasso)

bootstrap_pc_network<-function(samples_data_bootstrap,times_bootstrap,cut_off)
{
  adj_list = list()
  adj_final = NULL
  p_max_final = NULL
  p_max_list = list()
  for(i in 1:times_bootstrap)
  {
    new_sample_index = sample(seq(1,nrow(samples_data_bootstrap)), nrow(samples_data_bootstrap),replace = TRUE)
    new_samples = samples_data_bootstrap[new_sample_index,]
    # if(i==1) #or randomly choose
    # {
    #   save(new_samples, file=paste0(fold_path,"new_samples.Rdata"))
    # }
    corr_mat = pcor(new_samples) 
    suffStat <- list(C=corr_mat ,n=nrow(new_samples))
    alpha = 0.01   #alpha (p-value) threshold for pc function
    numCores = 16  #number of cores for stable.fast algo in pcalg
    m.max = 4      #max. number of conditioning variables
    pc_obj = pc(suffStat,
                indepTest=gaussCItest,
                p=ncol(corr_mat),
                skel.method="stable.fast",
                #fixedGaps = FixedGaps_mat,
                alpha = alpha,
                numCores = numCores ,
                m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209
                verbose= FALSE)
    # pc_obj = PC_Algorithm_GeneNetwork(new_samples)
    # ## estimate CPDAG and PDAG -- see  help(pc)
    # suffStat <- list(C = cor(new_samples), n = sample_num)
    # pc_obj <- pc(suffStat, indepTest = gaussCItest, p = node_num, alpha = 0.01)
    adj = as(pc_obj@graph,"matrix")
    adj_list[[i]] = adj
    
    p_max_list[[i]] = pc_obj@pMax
    
    if(i==1){adj_final = adj}
    else{ adj_final = adj_final + adj }
  }
  
  adj_final_output = adj_final
  library(popbio)
  p_max_final = mean.list(p_max_list) # initial value for p_max_final
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
    else
    { 
      adj_final_output[x_loc,y_loc] = 0 
      p_max_final[x_loc,y_loc] = -1
    }
  }
  graph_final_est = as(adj_final_output,"graphNEL")
  
  
  # for edge we final remain, calculate value for each of them
  non_zero_loc_edge_final = which(adj_final_output!=0,arr.ind = T)
  rows2 = dim(non_zero_loc_edge_final)[1] #num of row
  columns2 = dim(non_zero_loc_edge_final)[2] #num of column
  
  for (i in 1:rows2)
  {
    x_loc = non_zero_loc_edge_final[i,1]
    y_loc = non_zero_loc_edge_final[i,2]
    tmp_list_pvalue = list()
    tmp_count = 1 # count for tmp_list_pvalue
    for(j in 1:length(adj_list))
    {
      adj_list_each = adj_list[[j]]
      edge_mark = adj_list_each[x_loc,y_loc]
      if(edge_mark==1)
      {
        tmp_list_pvalue[[tmp_count]] = p_max_list[[j]][x_loc,y_loc]
        tmp_count = tmp_count+1
      }
    }
    p_max_final[x_loc,y_loc] = mean(unlist(tmp_list_pvalue)) 
  }
  
  bootstrap_result_list[[1]] = graph_final_est
  bootstrap_result_list[[2]] = p_max_final
  
  return(bootstrap_result_list)
}

bootstrap_pc_network2<-function(samples_data_bootstrap,times_bootstrap,cut_off)
{
  adj_list = list()
  adj_final = NULL
  p_max_final = NULL
  p_max_list = list()
  for(i in 1:times_bootstrap)
  {
    new_sample_index = sample(seq(1,nrow(samples_data_bootstrap)), nrow(samples_data_bootstrap),replace = TRUE)
    new_samples = samples_data_bootstrap[new_sample_index,]
    corr_mat = pcor(new_samples) 
    suffStat <- list(C=corr_mat ,n=nrow(new_samples))
    alpha = 0.01   #alpha (p-value) threshold for pc function
    numCores = 16  #number of cores for stable.fast algo in pcalg
    m.max = 4      #max. number of conditioning variables
    pc_obj = pc(suffStat,
                indepTest=gaussCItest,
                p=ncol(corr_mat),
                skel.method="stable.fast",
                #fixedGaps = FixedGaps_mat,
                alpha = alpha,
                numCores = numCores ,
                m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209
                verbose= FALSE)
    # pc_obj = PC_Algorithm_GeneNetwork(new_samples)
    # ## estimate CPDAG and PDAG -- see  help(pc)
    # suffStat <- list(C = cor(new_samples), n = sample_num)
    # pc_obj <- pc(suffStat, indepTest = gaussCItest, p = node_num, alpha = 0.01)
    adj = as(pc_obj@graph,"matrix")
    adj_list[[i]] = adj
    # p_max_list[[i]] = pc_obj@pMax
    
    adj_tmp = adj
    pMax_tmp = pc_obj@pMax
    pMax_new = pMax_tmp*adj_tmp
    pMax_new[adj==0]=-1
    p_max_list[[i]] = pMax_new
    if(i==1){
      adj_final = adj
    }
    else{ 
      adj_final = adj_final + adj 
      # p_max_final = p_max_final+pMax_new
      # p_max_list[[i]] = pc_obj@pMax
    }
  }
  
  adj_final_output = adj_final
  library(popbio)
  p_max_final = mean.list(p_max_list) # initial value for p_max_final
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
    else
    { 
      adj_final_output[x_loc,y_loc] = 0 
      p_max_final[x_loc,y_loc] = -1
      #cat("This is p_max_final",p_max_final[x_loc,y_loc],"\n")
    }
  }
  graph_final_est = as(adj_final_output,"graphNEL")
  
  bootstrap_result_list[[1]] = graph_final_est
  bootstrap_result_list[[2]] = p_max_final
  
  return(bootstrap_result_list)
}

bootstrap_pc_network3<-function(samples_data_bootstrap,times_bootstrap,cut_off)
{
  adj_list = list()
  adj_final = NULL # for adj=1
  adj_final2 = NULL # for adj=0
  # adj1 = NULL # for adj=1
  # adj2 = NULL # for adj=0
  p_max_final = NULL # for adj=1
  p_max_final2 = NULL # for adj=0
  p_max_list = list()
  
  for(i in 1:times_bootstrap)
  {
    new_sample_index = sample(seq(1,nrow(samples_data_bootstrap)), nrow(samples_data_bootstrap),replace = TRUE)
    new_samples = samples_data_bootstrap[new_sample_index,]
    # if(i==1) #or randomly choose
    # {
    #   save(new_samples, file=paste0(fold_path,"new_samples.Rdata"))
    # }
    # pc_obj = PC_Algorithm_GeneNetwork(new_samples)
    corr_mat = pcor(new_samples) 
    suffStat <- list(C=corr_mat ,n=nrow(new_samples))
    alpha = 0.01   #alpha (p-value) threshold for pc function
    numCores = 16  #number of cores for stable.fast algo in pcalg
    m.max = 4      #max. number of conditioning variables
    pc_obj = pc(suffStat,
                indepTest=gaussCItest,
                p=ncol(corr_mat),
                skel.method="stable.fast",
                #fixedGaps = FixedGaps_mat,
                alpha = alpha,
                numCores = numCores ,
                m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209
                verbose= FALSE)
    # ## estimate CPDAG and PDAG -- see  help(pc)
    # suffStat <- list(C = cor(new_samples), n = sample_num)
    # pc_obj <- pc(suffStat, indepTest = gaussCItest, p = node_num, alpha = 0.01)
    adj = as(pc_obj@graph,"matrix")
    
    adj_list[[i]] = adj
    p_max_list[[i]] = pc_obj@pMax
    pmax = pc_obj@pMax
    
    adj2 = adj
    adj2[adj2==1] = 100
    adj2[adj2==0] = 1
    adj2[adj2==100] = 0  # exchange 0 and 1 in adj
    
    if(i==1){
      adj_final = adj
      adj_final2 = adj2
      p_max_final = pmax*adj
      p_max_final2 = pmax*adj2
    }
    else{ 
      adj_final = adj_final + adj 
      adj_final2 = adj_final2 + adj2  
      p_max_final = p_max_final + pmax*adj
      p_max_final2 = p_max_final2 + pmax*adj2
    }
  }
  
  adj_final_output = adj_final
  library(popbio)
  # p_max_final = mean.list(p_max_list) # initial value for p_max_final
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
    else
    { 
      adj_final_output[x_loc,y_loc] = 0 
      p_max_final[x_loc,y_loc] = -1
    }
  }
  graph_final_est = as(adj_final_output,"graphNEL")
  
  p_max_output = adj_final_output 
  for(a in 1:nrow(adj_final_output))
  {
    for(b in 1:ncol(adj_final_output))
    {
      if(adj_final_output[a,b]==1)
      {
        p_max_output[a,b] = p_max_final[a,b]/adj_final[a,b]
      }
      else if(adj_final_output[a,b]==0)
      {
        p_max_output[a,b] = p_max_final2[a,b]/adj_final2[a,b]
      }
    }
  }
  
  bootstrap_result_list[[1]] = graph_final_est
  bootstrap_result_list[[2]] = p_max_output
  
  return(bootstrap_result_list)
}
#***************************************************************************************
# infer the gene-gene network with select genes by bootstrapped pc algorithm 
#***************************************************************************************
# The file directory for resid_mat with selected genes
resid_mat_selected_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Results_PCSelect_Based/"
#resid_mat_selected_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Result_Selected/"

files = list.files(path=resid_mat_selected_prefix)
# Iterative all files in the directory
for(i in 1:length(files)){
  #*************************************
  #  Step 1: load residual_mat 
  #*************************************
  filepath = paste0(resid_mat_selected_prefix,files[i])
  load(filepath) ## object: resid_mat_selected
  
  #*************************************
  #  Step 2: run bootstrapped pc 
  #*************************************
  bootstrap_result_list = list()
  pc_graph = bootstrap_pc_network3(resid_mat_selected,100,0.6)
  graph_est = pc_graph[[1]]
  pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_PCSelect/"
  file_name = gsub("_resid_mat_selected.Rdata","",files[i])
  save(pc_graph,file = paste0(pc_graph_prefix,file_name,"_Causal_Network_PCSelect.Rdata"))
}


