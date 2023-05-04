
# Written by yinly, May 11

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
library(Matrix)

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
  # Revised by yinly: August 26, 2020
  # Description: Calcuate E(V), i.e., expected number of falsely identified edges
  qlambda_sum = 0
  
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
    #***********************************************************************************
    # Revised by yinly, Oct 12, 2020
    # Description: update the pMax matrix derived from the pc algorithm accordingly
    #***********************************************************************************
    pMax_temp = pc_obj@pMax
    pMax_temp = forceSymmetric(pMax_temp,uplo="U")
    pMax_update = pMax_temp
    pMax_update[which(pMax_update==-1)] = 1
    pMax_final = as.matrix(pMax_update)
    pc_obj@pMax = pMax_final
    # end of revision
    
    adj = as(pc_obj@graph,"matrix")

    # Revised by yinly, August 26,2020
    # Description: update the qlambda_sum for each bootstrap iteration
    qlambda_sum = qlambda_sum + sum(adj)
    
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
  # Revised by yinly, August 26, 2020
  # Description: calculate E(V) 
  # Reference: https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2010.00740.x
  qlambda = qlambda_sum/times_bootstrap
  p = nrow(adj_final_output) * (nrow(adj_final_output)-1)
  EV = (1/(2*cut_off - 1)) * ((qlambda**2)/p) 

  bootstrap_result_list[[1]] = graph_final_est
  bootstrap_result_list[[2]] = p_max_output
  # Revised by yinly, August 26,2020
  # Description: include EV in the result list
  bootstrap_result_list[[3]] = EV
  bootstrap_result_list[[4]] = qlambda
  
  return(bootstrap_result_list)
}
#***************************************************************************************
# infer the gene-gene network with select genes by bootstrapped pc algorithm 
#***************************************************************************************
# The file directory for resid_mat with selected genes
#resid_mat_selected_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Results_PCSelect_Based/"
resid_mat_selected_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Residual_Result_Selected_CIS_Trans/"

files = list.files(path=resid_mat_selected_prefix)
EV_Result = matrix(0,nrow = length(files),ncol = 6)
colnames(EV_Result) = c("Trait","EV","qlamda","p","TypeI_error","FDR")
# Iterative all files in the directory
# Revised by yinly, 1 April 2022
# Description: perform the seletive analysis for covid patients with updated dataset
index <- c(55,56,89,90)
#for(i in 1:length(files)){
for(j in 1:length(index)){
  i <- index[j]
 #*************************************
 #  Step 1: load residual_mat 
 #*************************************
 #file = "Whole_Blood_CAD_resid_mat_selected.Rdata"
 filepath = paste0(resid_mat_selected_prefix,files[i])
 load(filepath) ## object: resid_mat_selected

 #*************************************
 #  Step 2: run bootstrapped pc 
 #*************************************
 bootstrap_result_list = list()
 pc_graph = bootstrap_pc_network3(resid_mat_selected,100,0.6)
 graph_est = pc_graph[[1]]
 # Revised by yinly, July 23, 2021
 # Description: calculate stability selection results for all traits
 pc_graph_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/GTEx_Selected_Update_Second/"
 file_name = gsub("_resid_mat_selected.Rdata","",files[i])
 #file_name = gsub("_resid_mat_selected.Rdata","",file)
 save(pc_graph,file = paste0(pc_graph_prefix,file_name,"_Causal_Network_PCSelect.Rdata"))
 # Revised by yinly: August 26, 2020
 # Description: output the EV of the inferred graph for all studied traits
 # Revised by yinly: July 23, 2021
 # Description: output EV, typeI error and FDR for all studied traits
 pc_graph_EV = pc_graph[[3]]
 pc_graph_qlambda <- pc_graph[[4]]
 pc_graph_p = nrow(pc_graph[[2]]) * (nrow(pc_graph[[2]]) - 1) 
 EV_Result[i,1] = file_name
 EV_Result[i,2] = pc_graph_EV
 EV_Result[i,3] <- pc_graph_qlambda
 EV_Result[i,4] = pc_graph_p
 EV_Result[i,5] <- pc_graph_EV/pc_graph_p
 EV_Result[i,6] <- pc_graph_EV/pc_graph_qlambda
 cat("The inference of GTEx Network for the",i,"trait",file_name,"is completed! \n")
}
EV_Result_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/"
fwrite(as.data.frame(EV_Result),paste0(EV_Result_prefix,"GTEx_Stability_Selection_Results_Update_COVID.csv"),sep="\t")
