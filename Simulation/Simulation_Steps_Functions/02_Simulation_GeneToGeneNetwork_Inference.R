#############################################################################
# Third Step : applied PC-algorithm on simulation network
# add bootstrap
#############################################################################
PC_Algorithm_GeneNetwork<-function(samples_data_bootstrap)
{
  
  resid_mat = samples_data_bootstrap
  
 
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
  
  cov_mat = coop::covar(resid_mat)
  
  t1=proc.time()
  
  rho_theoretical =  sqrt( log(ncol(resid_mat)) / nrow(resid_mat) )     ## https://arxiv.org/pdf/1403.6752.pdf p.1219 / SILGGM publication https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006369
  rho_theoretical = 1e-10
  print(rho_theoretical)
  library(glassoFast)
  
  glasso_obj = glassoFast(S=cov_mat,
                          rho= rho_theoretical)
  
  
  proc.time()-t1
  file_path = "/exeh_4/yaning_feng/04_Simulation/ForTest/"
  # save(glasso_obj, file=paste0(file_path,"Step2_glasso_file2.Rdata"))
  
  
  #####################Step 3: applied PC agrithomn#######################
  alpha = 0.01   #alpha (p-value) threshold for
  numCores = 16  #number of cores for stable.fast algo in pcalg
  m.max = 4      #max. number of conditioning variables
  
  library(CePa)
  library(ParallelPC)
  library(reshape2)
  library(data.table)
  library(dplyr)
  library(pcalg)
  library(coop)
  library(parallel)
  
  #*************************************
  #  pc algorithm applied to GTEx
  #***************************************
  corr_mat = pcor(resid_mat)
  #save(corr_mat, file="/exeh_4/sohc/network_GWAS/Correl_matrix_GTEx_RNA_seq_coronary_artery.Rdata")
  suffStat <- list(C=corr_mat ,n=nrow(resid_mat)) ##pcor is a function from "coop" for fast computation of correlation matrix
  
  #**********************
  # enforce symmetricity in the partial correlation matrix obtained from glasso
  #**********************************
  library(Matrix)
  FixedGaps_mat = forceSymmetric(glasso_obj$wi)
  FixedGaps_mat[FixedGaps_mat!=0] <- 999
  FixedGaps_mat[FixedGaps_mat==0] <- 1
  FixedGaps_mat[FixedGaps_mat==999] <- 0
  FixedGaps_mat = matrix( as.logical(FixedGaps_mat), nrow = nrow(FixedGaps_mat) )
  
  t1=proc.time()
  
  pc_obj = pc(suffStat,
              indepTest=gaussCItest,
              p=ncol(corr_mat),
              skel.method="stable.fast",
              # skel.method="stable",
              fixedGaps = FixedGaps_mat,
              alpha = alpha,
              numCores = numCores ,
              m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209
              verbose= FALSE)
  proc.time()-t1
  # There is an error in pMax when i>j, so we need adjustment.
  pmax = pc_obj@pMax
  adj = as(pc_obj@graph,"matrix")
  
  
  pmax_new = pmax
  # pmax_new[pmax_new==-1] = 1
  pmax_new = forceSymmetric(pmax_new)
  # pmax_new[which(adj==0 & pmax_new < alpha, arr.ind = T)] = 1
  
  # pc_obj@pMax = pmax_new
  pc_obj_list = list()
  pc_obj_list[[1]] = pc_obj@graph
  pc_obj_list[[2]] = pmax_new
  
  return(pc_obj_list)
}

PC_iterate<-function(cov_mat)
{
  t1=proc.time()
  
  rho_theoretical =  sqrt( log(ncol(resid_mat)) / nrow(resid_mat) )     ## https://arxiv.org/pdf/1403.6752.pdf p.1219 / SILGGM publication https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006369
  rho_theoretical = 1e-10
  print(rho_theoretical)
  library(glassoFast)
  
  glasso_obj = glassoFast(S=cov_mat,
                          rho= rho_theoretical)
  
  proc.time()-t1
  file_path = "/exeh_4/yaning_feng/04_Simulation/ForTest/"
  
  
  #####################Step 3: applied PC agrithomn#######################
  alpha = alpha  #alpha (p-value) threshold for
  numCores = 16  #number of cores for stable.fast algo in pcalg
  m.max = 4      #max. number of conditioning variables
  
  library(CePa)
  library(ParallelPC)
  library(reshape2)
  library(data.table)
  library(dplyr)
  library(pcalg)
  library(coop)
  library(parallel)
  
  #*************************************
  #  pc algorithm applied to GTEx
  #***************************************

  corr_mat = cov2cor(cov_mat)
  #save(corr_mat, file="/exeh_4/sohc/network_GWAS/Correl_matrix_GTEx_RNA_seq_coronary_artery.Rdata")
  suffStat <- list(C=corr_mat ,n=nrow(resid_mat)) ##pcor is a function from "coop" for fast computation of correlation matrix
  
  #**********************
  # enforce symmetricity in the partial correlation matrix obtained from glasso
  #**********************************
  library(Matrix)
  FixedGaps_mat = forceSymmetric(glasso_obj$wi)
  FixedGaps_mat[FixedGaps_mat!=0] <- 999
  FixedGaps_mat[FixedGaps_mat==0] <- 1
  FixedGaps_mat[FixedGaps_mat==999] <- 0
  FixedGaps_mat = matrix( as.logical(FixedGaps_mat), nrow = nrow(FixedGaps_mat) )
  
  t1=proc.time()
  
  pc_obj = pc(suffStat,
              indepTest=gaussCItest,
              p=ncol(corr_mat),
              skel.method="stable.fast",
              fixedGaps = FixedGaps_mat,
              alpha = 0.01,
              numCores = numCores ,
              m.max = m.max ,  #ref: reduced PC paper https://arxiv.org/abs/1806.06209
              verbose= FALSE)
  proc.time()-t1

  return(pc_obj)
}

bootstrap_pc_network3<-function(samples_data_bootstrap,times_bootstrap,cut_off)
{
  adj_list = list()
  adj_final = NULL # for adj=1
  adj_final2 = NULL # for adj=0
  
  p_max_final = NULL # for adj=1
  p_max_final2 = NULL # for adj=0
  p_max_list = list()
  
  for(i in 1:times_bootstrap)
  {
    new_sample_index = sample(seq(1,nrow(samples_data_bootstrap)), nrow(samples_data_bootstrap),replace = TRUE)
    new_samples = samples_data_bootstrap[new_sample_index,]
   
    pc_obj_list = PC_Algorithm_GeneNetwork(new_samples)
    
    adj = as(pc_obj_list[[1]],"matrix")
    
    adj_list[[i]] = adj
    p_max_list[[i]] = pc_obj_list[[2]]
    pmax = pc_obj_list[[2]]
    
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
  bootstrap_result_list = list()
  bootstrap_result_list[[1]] = graph_final_est
  bootstrap_result_list[[2]] = p_max_output
  
  return(bootstrap_result_list)
}



