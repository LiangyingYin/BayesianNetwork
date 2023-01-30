source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Extension_Functions/Simulation_ExtensionFunction_Mediator_functions.R")

GetMoreDetail_OverallPerformance<-function(dt_IDA)
{
  #obtain minimum jointIDA, result in est_jointIDA
  est_jointIDA = NULL
  ind1 = which(colnames(dt_IDA)=="true_jointIDA") # obtain the column No of true_jointIDA
  ind2 = which(colnames(dt_IDA)=="Rank_ABS_true_IDA") # obtain the column No of Rank_ABS_true_IDA
  jointIDA_cols = as.data.frame( dt_IDA[,(ind1+1):(ind2-1)], nrow= nrow(dt_IDA)) # all the alternative value of tureIDA
  jointIDA_cols[jointIDA_cols=="NA"] <- 99999
  for (j in 1:nrow(jointIDA_cols) ) {
    rowi = as.numeric(as.character(jointIDA_cols[j,]))
    min_index = which.min(abs(rowi))
    est_jointIDA[j] = jointIDA_cols[j,min_index]
  }
  est_jointIDA = as.numeric(as.character(est_jointIDA))
  dt_IDA = data.frame(dt_IDA, est_jointIDA)
  
  
  #obtain minimum IDA, result in Estimate_IDA_final
  Estimate_IDA_final = NULL
  ind1 = which(colnames(dt_IDA)=="Estimate_IDA")
  ind2 = which(colnames(dt_IDA)=="true_jointIDA")
  IDA_cols = as.data.frame( dt_IDA[,ind1:(ind2-1)], nrow= nrow(dt_IDA)) 
  IDA_cols[IDA_cols=="NA"] <- 99999
  for (j in 1:nrow(IDA_cols) ) {
    rowi = as.numeric(as.character(IDA_cols[j,]))
    min_index = which.min(abs(rowi))
    Estimate_IDA_final[j] = IDA_cols[j,min_index]
  }
  Estimate_IDA_final = as.numeric(as.character(Estimate_IDA_final))	
  dt_IDA = data.frame(dt_IDA, Estimate_IDA_final)
  return(dt_IDA)
}

Get_related_variable<-function(dt_IDA,adj,g1,g2,g3,gene_index)
{
  
  #save correlation of IDA
  Cor_TureIDA_EstIDA_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final,  use='complete.obs' ) 
  Cor_TureIDA_EstIDA_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final, use='complete.obs', method="spearman")
  Cor_TureIDA_ImputeGene_Rcof_Pearson	= cor(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof, use='complete.obs')
  Cor_TureIDA_ImputeGene_Rcof_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof, use='complete.obs', method="spearman")
  Cor_TureIDA_rowGene_Rcof_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof)
  Cor_TureIDA_rowGene_Rcof_Spearman = cor(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof, method="spearman")
  
  
  #save correlation of nonzero IDA
  ind_IDA =  which(abs(dt_IDA$true_IDA)>1e-10)
  Cor_TureIDA_EstIDA_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Estimate_IDA_final[ind_IDA],  use='complete.obs' )
  Cor_TureIDA_EstIDA_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Estimate_IDA_final[ind_IDA], use='complete.obs', method="spearman")
  Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_imputedata_regression_Cof[ind_IDA], use='complete.obs' )
  Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_imputedata_regression_Cof[ind_IDA], use='complete.obs', method="spearman" )
  Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_rowdata_regression_Cof[ind_IDA])
  Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_IDA[ind_IDA], dt_IDA$Gene_rowdata_regression_Cof[ind_IDA], method="spearman")
  
  #save RMSE
  RMSE_TrueIDA_EstIDA = RMSE(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final)
  RMSE_TrueIDA_ImputeGene_Rcof = RMSE(dt_IDA$true_IDA, dt_IDA$Gene_imputedata_regression_Cof)
  RMSE_TrueIDA_rowGene_Rcof = RMSE(dt_IDA$true_IDA, dt_IDA$Gene_rowdata_regression_Cof)
  
  
  #____________________________
  # correlation of joint IDA
  #____________________________
  dt_IDA$true_jointIDA <- as.numeric(as.character(dt_IDA$true_jointIDA))
  # save correlation of joint IDA
  Cor_TureJointIDA_EstJointIDA_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  Cor_TureJointIDA_EstJointIDA_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA, method="spearman")
  
  # modified by fyn 14/2/2020
  Cor_TureJointIDA_ImputeGene_Rcof_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof)
  Cor_TureJointIDA_ImputeGene_Rcof_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof,  method="spearman")
  
  Cor_TureJointIDA_rowGene_Rcof_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof)
  Cor_TureJointIDA_rowGene_Rcof_Spearman = cor(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof,  method="spearman")
  
  
  #________________________________________________
  # correlation of joint IDA which are non-zero
  #________________________________________________
  dt_IDA$true_jointIDA <- as.numeric(as.character(dt_IDA$true_jointIDA))
  ind_nonzero =  which(abs(as.numeric(dt_IDA$true_jointIDA))>1e-10)
  
  Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$est_jointIDA[ind_nonzero])
  Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$est_jointIDA[ind_nonzero], method="spearman")
  
  Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_imputedata_regression_Cof[ind_nonzero])
  Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_imputedata_regression_Cof[ind_nonzero],method="spearman")
  
  Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_rowdata_regression_Cof[ind_nonzero])	
  Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero = cor(dt_IDA$true_jointIDA[ind_nonzero], dt_IDA$Gene_rowdata_regression_Cof[ind_nonzero],  method="spearman")
  
  RMSE_TrueJointIDA_EstJointIDA = RMSE(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  
  RMSE_TrueJointIDA_ImputeGene_Rcof = RMSE(dt_IDA$true_jointIDA, dt_IDA$Gene_imputedata_regression_Cof)
  RMSE_TrueJointIDA_rowGene_Rcof = RMSE(dt_IDA$true_jointIDA, dt_IDA$Gene_rowdata_regression_Cof)
  
  # proportion of gene-gene connection
  gene_mat = adj[gene_index,gene_index]
  sum(gene_mat!=0)/nrow(gene_mat)^2
  
  proportion_GeneToGene_Net = sum(gene_mat!=0)/nrow(gene_mat)^2
  
  # correlation of regression coef derived from real expr vs imputed expr
  Cor_rowGene_ImputeGene_Rcof = cor(dt_IDA$Gene_rowdata_regression_Cof, dt_IDA$Gene_imputedata_regression_Cof)
  
  shd_ture_estGrapph = shd(g1,g2)
  shd_ture_randomGraph = shd(g1,g3)
  
  # correlation and MSE
  related_variable = c(Cor_TureIDA_EstIDA_Pearson,Cor_TureIDA_EstIDA_Spearman,Cor_TureIDA_ImputeGene_Rcof_Pearson,Cor_TureIDA_ImputeGene_Rcof_Spearman,	Cor_TureIDA_rowGene_Rcof_Pearson,Cor_TureIDA_rowGene_Rcof_Spearman,Cor_TureIDA_EstIDA_Pearson_Nonzero,	Cor_TureIDA_EstIDA_Spearman_Nonzero,	Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero,Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero,	Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero,	Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero,	RMSE_TrueIDA_EstIDA,	RMSE_TrueIDA_ImputeGene_Rcof,	RMSE_TrueIDA_rowGene_Rcof,	Cor_TureJointIDA_EstJointIDA_Pearson,	Cor_TureJointIDA_EstJointIDA_Spearman,	Cor_TureJointIDA_ImputeGene_Rcof_Pearson,	Cor_TureJointIDA_ImputeGene_Rcof_Spearman,	Cor_TureJointIDA_rowGene_Rcof_Pearson,	Cor_TureJointIDA_rowGene_Rcof_Spearman,	Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero,	Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero,	Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero,	Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero,	Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero,	Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero,	RMSE_TrueJointIDA_EstJointIDA ,	RMSE_TrueJointIDA_ImputeGene_Rcof ,	RMSE_TrueJointIDA_rowGene_Rcof ,	proportion_GeneToGene_Net,	Cor_rowGene_ImputeGene_Rcof,	shd_ture_estGrapph,	shd_ture_randomGraph)
  related_variable_colname = c("Cor_TureIDA_EstIDA_Pearson","Cor_TureIDA_EstIDA_Spearman","Cor_TureIDA_ImputeGene_Rcof_Pearson","Cor_TureIDA_ImputeGene_Rcof_Spearman","Cor_TureIDA_rowGene_Rcof_Pearson","Cor_TureIDA_rowGene_Rcof_Spearman","Cor_TureIDA_EstIDA_Pearson_Nonzero","Cor_TureIDA_EstIDA_Spearman_Nonzero","Cor_TureIDA_ImputeGene_Rcof_Pearson_Nonzero","Cor_TureIDA_ImputeGene_Rcof_Spearman_Nonzero","Cor_TureIDA_rowGene_Rcof_Pearson_Nonzero","Cor_TureIDA_rowGene_Rcof_Spearman_Nonzero","RMSE_TrueIDA_EstIDA","RMSE_TrueIDA_ImputeGene_Rcof","RMSE_TrueIDA_rowGene_Rcof","Cor_TureJointIDA_EstJointIDA_Pearson","Cor_TureJointIDA_EstJointIDA_Spearman","Cor_TureJointIDA_ImputeGene_Rcof_Pearson","Cor_TureJointIDA_ImputeGene_Rcof_Spearman","Cor_TureJointIDA_rowGene_Rcof_Pearson","Cor_TureJointIDA_rowGene_Rcof_Spearman","Cor_TureJointIDA_EstJointIDA_Pearson_Nonzero","Cor_TureJointIDA_EstJointIDA_Spearman_Nonzero","Cor_TureJointIDA_ImputeGene_Rcof_Pearson_Nonzero","Cor_TureJointIDA_ImputeGene_Rcof_Spearman_Nonzero","Cor_TureJointIDA_rowGene_Rcof_Pearson_Nonzero","Cor_TureJointIDA_rowGene_Rcof_Spearman_Nonzero","RMSE_TrueJointIDA_EstJointIDA ","RMSE_TrueJointIDA_ImputeGene_Rcof","RMSE_TrueJointIDA_rowGene_Rcof","proportion_GeneToGene_Net","Cor_rowGene_ImputeGene_Rcof","shd_ture_estGrapph","shd_ture_randomGraph")
  list_related_variable = list()
  list_related_variable[["related_variable"]] = related_variable
  list_related_variable[["related_variable_colname"]] = related_variable_colname
  
  return(list_related_variable)
}

Performance_mediatorDetect_FromGraphStructure<-function(ifmediator_gene_true,ifmediator_gene_est,genes_num)
{
  ######################################################################################
  # Find the mediator in the gene-outcome graph by using graph structure.
  ifmediator_gene_true_std = Standardize_Mediator_result(ifmediator_gene_true)
  ifmediator_gene_est_std = Standardize_Mediator_result(ifmediator_gene_est)
  
  ifmediator_gene_true = ifmediator_gene_true[1:genes_num]
  ifmediator_gene_est = ifmediator_gene_est[1:genes_num]
  
  tab <- table(ifmediator_gene_est_std,ifmediator_gene_true_std)[2:1,2:1]
  mediator_rval <- epi.tests(tab, conf.level = 0.95)
  mediator_rval_aprev = mediator_rval$rval['aprev'][[1]][1,1] # apparent prevalence.
  mediator_rval_tprev = mediator_rval$rval['tprev'][[1]][1,1]  # true prevalence.
  mediator_rval_se = mediator_rval$rval['se'][[1]][1,1]   # Sensitivity
  mediator_rval_sp = mediator_rval$rval['sp'][[1]][1,1]   # Specificity
  mediator_rval_ppv = mediator_rval$rval['ppv'][[1]][1,1]   # Positive predictive value
  mediator_rval_npv = mediator_rval$rval['npv'][[1]][1,1]   # Negative predictive value
  mediator_rval_plr = mediator_rval$rval['plr'][[1]][1,1]   # Positive likelihood ratio
  mediator_rval_nlr = mediator_rval$rval['nlr'][[1]][1,1]   # Negative likelihood ratio
  
  mediator_rval_combine = c(mediator_rval_aprev,mediator_rval_tprev,mediator_rval_se,mediator_rval_sp,mediator_rval_ppv,mediator_rval_npv,mediator_rval_plr,mediator_rval_nlr)
  mediator_rval_combine_colname = c("mediator_rval_aprev","mediator_rval_tprev","mediator_rval_se","mediator_rval_sp","mediator_rval_ppv","mediator_rval_npv","mediator_rval_plr","mediator_rval_nlr")
  
  list_mediator_rval_combine = list()
  list_mediator_rval_combine[["mediator_rval_combine"]] = mediator_rval_combine
  list_mediator_rval_combine[["mediator_rval_combine_colname"]] = mediator_rval_combine_colname
  
  return(list_mediator_rval_combine)
}

PC_Simple_Performance<-function(adj,gene_index,n3,pcSimple.fi)
{
  #_________________________________________________
  # this part examines the performance of PC simple in detecting direct causal variants
  #_____________________________________________________
  Pheno_vec = adj[gene_index,n3]
  
  gene_list = attributes(pcSimple.fit$G)
  pcsimple_index = attributes (which(pcSimple.fit$G==TRUE))$names # fyn: gene index with true relations to phenotype
  pcsimple_index = as.numeric(pcsimple_index) 
  real_directCausal = gene_index[Pheno_vec!=0]
  
  cell_11 = length(intersect(pcsimple_index, real_directCausal))
  cell_12 = length(pcsimple_index) - cell_11
  cell_21 = length(real_directCausal) - cell_11
  cell_22 = MyData[x,2] - cell_11 - cell_12 - cell_21  # MyData[x,2] is gene_num
  
  PCsimple_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  
  
  return(PCsimple_mat)
}

Impute_Performance<-function(dt_IDA,MyData,Pheno_vec)
{
  # this part examines the performance of PrediXcan (TWAS) in detecting direct causal variants
  #_____________________________________________________
  impExpr_ind = which( dt_IDA$Gene_imputedata_regression_pval <= 0.001)
  real_ind = which(Pheno_vec!=0) 
  cell_11 = length(  intersect(impExpr_ind, real_ind)  )
  cell_12 = length(impExpr_ind) - cell_11
  cell_21 = length(real_ind ) - cell_11
 
  cell_22 = MyData[x,2] - cell_11 - cell_12 - cell_21
  ImpExpr_mat = matrix( c(cell_11, cell_12, cell_21, cell_22), byrow=T, nrow=2) 
  return(ImpExpr_mat)
}

Get_epi_fisher_Summary<-function(PCsimple_mat,ImpExpr_mat)
{
  epi.tests_PCsimple = epi.tests(PCsimple_mat)
  fisher.test_PCsimple = fisher.test(PCsimple_mat)
  
  epi.tests_ImpExpr = epi.tests(ImpExpr_mat)
  fisher.test_ImpExpr = fisher.test(ImpExpr_mat)
  
  epi_PCsimple_aprev = epi.tests_PCsimple$rval['aprev'][[1]] # apparent prevalence.
  epi_PCsimple_tprev = epi.tests_PCsimple$rval['tprev'][[1]]  # true prevalence.
  epi_PCsimple_se = epi.tests_PCsimple$rval['se'][[1]]   # Sensitivity
  epi_PCsimple_sp = epi.tests_PCsimple$rval['sp'][[1]]   # Specificity
  epi_PCsimple_ppv = epi.tests_PCsimple$rval['ppv'][[1]]   # Positive predictive value
  epi_PCsimple_npv = epi.tests_PCsimple$rval['npv'][[1]]   # Negative predictive value
  epi_PCsimple_plr = epi.tests_PCsimple$rval['plr'][[1]]   # Positive likelihood ratio
  epi_PCsimple_nlr = epi.tests_PCsimple$rval['nlr'][[1]]   # Negative likelihood ratio
  
  epi_ImpExpr_aprev = epi.tests_ImpExpr$rval['aprev'][[1]] # apparent prevalence.
  epi_ImpExpr_tprev = epi.tests_ImpExpr$rval['tprev'][[1]]  # true prevalence.
  epi_ImpExpr_se = epi.tests_ImpExpr$rval['se'][[1]]   # Sensitivity
  epi_ImpExpr_sp = epi.tests_ImpExpr$rval['sp'][[1]]   # Specificity
  epi_ImpExpr_ppv = epi.tests_ImpExpr$rval['ppv'][[1]]   # Positive predictive value
  epi_ImpExpr_npv = epi.tests_ImpExpr$rval['npv'][[1]]   # Negative predictive value
  epi_ImpExpr_plr = epi.tests_ImpExpr$rval['plr'][[1]]   # Positive likelihood ratio
  epi_ImpExpr_nlr = epi.tests_ImpExpr$rval['nlr'][[1]]   # Negative likelihood ratio
  
  ##extract each data in epi.test
  epi_PCsimple_aprev_est = epi_PCsimple_aprev[1,1]
  epi_PCsimple_aprev_lower = epi_PCsimple_aprev[1,2]
  epi_PCsimple_aprev_upper = epi_PCsimple_aprev[1,3]
  epi_PCsimple_tprev_est = epi_PCsimple_tprev[1,1]
  epi_PCsimple_tprev_lower = epi_PCsimple_tprev[1,2]
  epi_PCsimple_tprev_upper = epi_PCsimple_tprev[1,3]
  epi_PCsimple_se_est = epi_PCsimple_se[1,1]
  epi_PCsimple_se_lower = epi_PCsimple_se[1,2]
  epi_PCsimple_se_upper = epi_PCsimple_se[1,3]
  epi_PCsimple_sp_est = epi_PCsimple_sp[1,1]
  epi_PCsimple_sp_lower = epi_PCsimple_sp[1,2]
  epi_PCsimple_sp_upper = epi_PCsimple_sp[1,3]
  epi_PCsimple_ppv_est = epi_PCsimple_ppv[1,1]
  epi_PCsimple_ppv_lower = epi_PCsimple_ppv[1,2]
  epi_PCsimple_ppv_upper = epi_PCsimple_ppv[1,3]
  epi_PCsimple_npv_est = epi_PCsimple_npv[1,1]
  epi_PCsimple_npv_lower = epi_PCsimple_npv[1,2]
  epi_PCsimple_npv_upper = epi_PCsimple_npv[1,3]
  epi_PCsimple_plr_est = epi_PCsimple_plr[1,1]
  epi_PCsimple_plr_lower = epi_PCsimple_plr[1,2]
  epi_PCsimple_plr_upper = epi_PCsimple_plr[1,3]
  epi_PCsimple_nlr_est = epi_PCsimple_nlr[1,1]
  epi_PCsimple_nlr_lower = epi_PCsimple_nlr[1,2]
  epi_PCsimple_nlr_upper = epi_PCsimple_nlr[1,3]
  
  epi_ImpExpr_aprev_est = epi_ImpExpr_aprev[1,1]
  epi_ImpExpr_aprev_lower = epi_ImpExpr_aprev[1,2]
  epi_ImpExpr_aprev_upper = epi_ImpExpr_aprev[1,3]
  epi_ImpExpr_tprev_est = epi_ImpExpr_tprev[1,1]
  epi_ImpExpr_tprev_lower = epi_ImpExpr_tprev[1,2]
  epi_ImpExpr_tprev_upper = epi_ImpExpr_tprev[1,3]
  epi_ImpExpr_se_est = epi_ImpExpr_se[1,1]
  epi_ImpExpr_se_lower = epi_ImpExpr_se[1,2]
  epi_ImpExpr_se_upper = epi_ImpExpr_se[1,3]
  epi_ImpExpr_sp_est = epi_ImpExpr_sp[1,1]
  epi_ImpExpr_sp_lower = epi_ImpExpr_sp[1,2]
  epi_ImpExpr_sp_upper = epi_ImpExpr_sp[1,3]
  epi_ImpExpr_ppv_est = epi_ImpExpr_ppv[1,1]
  epi_ImpExpr_ppv_lower = epi_ImpExpr_ppv[1,2]
  epi_ImpExpr_ppv_upper = epi_ImpExpr_ppv[1,3]
  epi_ImpExpr_npv_est = epi_ImpExpr_npv[1,1]
  epi_ImpExpr_npv_lower = epi_ImpExpr_npv[1,2]
  epi_ImpExpr_npv_upper = epi_ImpExpr_npv[1,3]
  epi_ImpExpr_plr_est = epi_ImpExpr_plr[1,1]
  epi_ImpExpr_plr_lower = epi_ImpExpr_plr[1,2]
  epi_ImpExpr_plr_upper = epi_ImpExpr_plr[1,3]
  epi_ImpExpr_nlr_est = epi_ImpExpr_nlr[1,1]
  epi_ImpExpr_nlr_lower = epi_ImpExpr_nlr[1,2]
  epi_ImpExpr_nlr_upper = epi_ImpExpr_nlr[1,3]
  
  
  fisher_PCsimple_pvalue = fisher.test(PCsimple_mat)$p.value
  fisher_PCsimple_95PCI_1 = fisher.test(PCsimple_mat)$conf.int[1] # 95 percent confidence interval
  fisher_PCsimple_95PCI_2 = fisher.test(PCsimple_mat)$conf.int[2] # 95 percent confidence interval
  fisher_PCsimple_OR = as.numeric(fisher.test(PCsimple_mat)$estimate) # odds ratio
  
  fisher_ImpExpr_pvalue = fisher.test(ImpExpr_mat)$p.value
  fisher_ImpExpr_95PCI_1 = fisher.test(ImpExpr_mat)$conf.int[1] # 95 percent confidence interval
  fisher_ImpExpr_95PCI_2 = fisher.test(ImpExpr_mat)$conf.int[2] # 95 percent confidence interval
  fisher_ImpExpr_OR = as.numeric(fisher.test(ImpExpr_mat)$estimate) # odds ratio
  
  epi_fisher_testresult_colname = c("epi_PCsimple_aprev_est ","epi_PCsimple_tprev_est ","epi_PCsimple_se_est ","epi_PCsimple_sp_est ","epi_PCsimple_ppv_est ","epi_PCsimple_npv_est ","epi_PCsimple_plr_est ","epi_PCsimple_nlr_est ","epi_ImpExpr_aprev_est ","epi_ImpExpr_tprev_est ","epi_ImpExpr_se_est ","epi_ImpExpr_sp_est ","epi_ImpExpr_ppv_est ","epi_ImpExpr_npv_est ","epi_ImpExpr_plr_est ","epi_ImpExpr_nlr_est","fisher_PCsimple_pvalue ","fisher_PCsimple_95PCI_1 ","fisher_PCsimple_95PCI_2 ","fisher_PCsimple_OR","fisher_ImpExpr_pvalue","fisher_ImpExpr_95PCI_1","fisher_ImpExpr_95PCI_2","fisher_ImpExpr_OR","epi_PCsimple_aprev_lower ","epi_PCsimple_aprev_upper","epi_PCsimple_tprev_lower ","epi_PCsimple_tprev_upper","epi_PCsimple_se_lower ","epi_PCsimple_se_upper ","epi_PCsimple_sp_lower ","epi_PCsimple_sp_upper ","epi_PCsimple_ppv_lower ","epi_PCsimple_ppv_upper ","epi_PCsimple_npv_lower ","epi_PCsimple_npv_upper ","epi_PCsimple_plr_lower ","epi_PCsimple_plr_upper","epi_PCsimple_nlr_lower ","epi_PCsimple_nlr_upper ","epi_ImpExpr_aprev_lower","epi_ImpExpr_aprev_upper","epi_ImpExpr_tprev_lower","epi_ImpExpr_tprev_upper","epi_ImpExpr_se_lower","epi_ImpExpr_se_upper","epi_ImpExpr_sp_lower","epi_ImpExpr_sp_upper","epi_ImpExpr_ppv_lower","epi_ImpExpr_ppv_upper","epi_ImpExpr_npv_lower","epi_ImpExpr_npv_upper","epi_ImpExpr_plr_lower","epi_ImpExpr_plr_upper","epi_ImpExpr_nlr_lower","epi_ImpExpr_nlr_upper")
  epi_fisher_testresult = c(epi_PCsimple_aprev_est ,	epi_PCsimple_tprev_est ,	epi_PCsimple_se_est ,	epi_PCsimple_sp_est ,	epi_PCsimple_ppv_est ,	epi_PCsimple_npv_est ,	epi_PCsimple_plr_est ,	epi_PCsimple_nlr_est ,	epi_ImpExpr_aprev_est ,	epi_ImpExpr_tprev_est ,	epi_ImpExpr_se_est ,	epi_ImpExpr_sp_est ,	epi_ImpExpr_ppv_est ,	epi_ImpExpr_npv_est ,	epi_ImpExpr_plr_est ,	epi_ImpExpr_nlr_est,	fisher_PCsimple_pvalue ,	fisher_PCsimple_95PCI_1 ,	fisher_PCsimple_95PCI_2 ,	fisher_PCsimple_OR,	fisher_ImpExpr_pvalue,	fisher_ImpExpr_95PCI_1,	fisher_ImpExpr_95PCI_2,	fisher_ImpExpr_OR,	epi_PCsimple_aprev_lower ,	epi_PCsimple_aprev_upper,	epi_PCsimple_tprev_lower ,	epi_PCsimple_tprev_upper,	epi_PCsimple_se_lower ,	epi_PCsimple_se_upper ,	epi_PCsimple_sp_lower ,	epi_PCsimple_sp_upper ,	epi_PCsimple_ppv_lower ,	epi_PCsimple_ppv_upper ,	epi_PCsimple_npv_lower ,	epi_PCsimple_npv_upper ,	epi_PCsimple_plr_lower ,	epi_PCsimple_plr_upper,	epi_PCsimple_nlr_lower ,	epi_PCsimple_nlr_upper ,	epi_ImpExpr_aprev_lower,	epi_ImpExpr_aprev_upper,	epi_ImpExpr_tprev_lower,	epi_ImpExpr_tprev_upper,	epi_ImpExpr_se_lower,	epi_ImpExpr_se_upper,	epi_ImpExpr_sp_lower,	epi_ImpExpr_sp_upper,	epi_ImpExpr_ppv_lower,	epi_ImpExpr_ppv_upper,	epi_ImpExpr_npv_lower,	epi_ImpExpr_npv_upper,	epi_ImpExpr_plr_lower,	epi_ImpExpr_plr_upper,	epi_ImpExpr_nlr_lower,	epi_ImpExpr_nlr_upper)
  
  list_epi_fisher_testresult = list()
  list_epi_fisher_testresult[["epi_fisher_testresult_colname"]] = epi_fisher_testresult_colname
  list_epi_fisher_testresult[["epi_fisher_testresult"]] = epi_fisher_testresult
  
  return(list_epi_fisher_testresult)
}

Get_related_variable_pvalue_opt<-function(dt_IDA,dt_IDA2,dt_IDA5)
{
  
  #save correlation above
  Cor_TureIDA_EstIDA_Pearson = cor(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final,  use='complete.obs' ) 
  Cor_TureJointIDA_EstJointIDA_Pearson = cor(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  RMSE_TrueIDA_EstIDA = RMSE(dt_IDA$true_IDA, dt_IDA$Estimate_IDA_final)
  RMSE_TrueJointIDA_EstJointIDA = RMSE(dt_IDA$true_jointIDA, dt_IDA$est_jointIDA)
  
  #add info from pvalue optimizaiton for cov
  Vx2 = grep("true_jointIDA", colnames(dt_IDA2))+1
  Cor_TureIDA_EstIDA_Pearson_pvalue  = cor(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
  RMSE_TrueIDA_EstIDA_pvalue = RMSE(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
  Cor_TureJointIDA_EstJointIDA_Pearson_pvalue = cor(as.numeric(dt_IDA2$true_jointIDA),dt_IDA2[,Vx2])
  RMSE_TrueJointIDA_EstJointIDA_pvalue = RMSE(as.numeric(dt_IDA2$true_jointIDA),dt_IDA2[,Vx2])
  
  
  Vx2 = grep("true_jointIDA", colnames(dt_IDA5))+1
  Cor_TureIDA_EstIDA_Pearson_pvalue_iterate  = cor(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
  RMSE_TrueIDA_EstIDA_pvalue_iterate = RMSE(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
  Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate = cor(as.numeric(dt_IDA5$true_jointIDA),dt_IDA5[,Vx2])
  RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate = RMSE(as.numeric(dt_IDA5$true_jointIDA),dt_IDA5[,Vx2])
  
  related_variable_pvalue_opt = c(Cor_TureIDA_EstIDA_Pearson,Cor_TureIDA_EstIDA_Pearson_pvalue,Cor_TureIDA_EstIDA_Pearson_pvalue_iterate,Cor_TureJointIDA_EstJointIDA_Pearson,Cor_TureJointIDA_EstJointIDA_Pearson_pvalue,Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate ,RMSE_TrueIDA_EstIDA,RMSE_TrueIDA_EstIDA_pvalue,RMSE_TrueIDA_EstIDA_pvalue_iterate,RMSE_TrueJointIDA_EstJointIDA,RMSE_TrueJointIDA_EstJointIDA_pvalue,RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate)
  related_variable_pvalue_opt_colname = c("Cor_TureIDA_EstIDA_Pearson","Cor_TureIDA_EstIDA_Pearson_pvalue","Cor_TureIDA_EstIDA_Pearson_pvalue_iterate","Cor_TureJointIDA_EstJointIDA_Pearson","Cor_TureJointIDA_EstJointIDA_Pearson_pvalue","Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate","RMSE_TrueIDA_EstIDA","RMSE_TrueIDA_EstIDA_pvalue","RMSE_TrueIDA_EstIDA_pvalue_iterate","RMSE_TrueJointIDA_EstJointIDA","RMSE_TrueJointIDA_EstJointIDA_pvalue","RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate")
  
  list_related_variable_pvalue_opt = list()
  list_related_variable_pvalue_opt[["related_variable_pvalue_opt"]] = related_variable_pvalue_opt
  list_related_variable_pvalue_opt[["related_variable_pvalue_opt_colname"]] = related_variable_pvalue_opt_colname
  
  return(list_related_variable_pvalue_opt)
}
