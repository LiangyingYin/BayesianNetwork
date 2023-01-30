######################################################################################
# Load the source code and function code
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Extension_Functions/Simulation_ExtensionFunction_Fisher_IDA_topK.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Extension_Functions/Simulation_ExtensionFunction_Bootstrap.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Extension_Functions/Simulation_ExtensionFunction_Mediator_functions.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Extension_Functions/Simulation_ExtensionFunction_pcSelect_stableVer.R")


source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/01_Simulation_GenerateSamples.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/02_Simulation_GeneToGeneNetwork_Inference.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/03_Simulation_GeneToPhenoNetwork_Inference.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/04_Simulation_Merge.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/05_Simulation_UnivariateTest.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/06_Simulation_StrategyForAdjustment.R")
source("/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/Simulation_Steps_Functions/07_Simulation_OrganizeSummarizeResults.R")

######################################################################################
# Import the library
library(pcalg)
library(CePa)
library(reshape2)
library(data.table)
library(dplyr)
library(coop)
library(parallel)
library(graph)
library(grid)
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
library(matrixcalc)
library(CVglasso)
library(iterators)
library(foreach)
library(doParallel)
library(rlist)

######################################################################################
#Setting the source &results path and some Parameters
csv_path =  "/exeh_4/yaning_feng/04_Simulation/FinalVersion/Rundata/Re-OrganizeProgram/"
csv_filename = "Simulation_Overall_Rundata.csv"
result_path = csv_path

min_weight_gene_to_pheno = 0
max_weight_gene_to_pheno = 3
max_weight_snp_to_gene = 0.5 #in real situation, relations between snps and genes is weak  
prop_snp_related_gene = 0.01  #add by fyn 27/11/2020. set the proportion of snps related to gene 
prop_snp_predict_gene = 1

cv.glmnet_nfolds = 5  # for cross-validation
eMat_g_tmp_rnorm_sd = 2  #gene expression data error variance
prop_direct_causal_genes = 1/5

#some parameters about bootstrap
times_bootstrap = 100
cut_off = 0.6


######################################################################################
#Getting the scenario info.
MyData <- read.csv(file=paste0(csv_path,csv_filename), header=TRUE, sep=",")
rows_testdata = dim(MyData)[1]


for (x in 1:rows_testdata)
{
  ######################################################################################
  #Create the num variable and list.
  list_dt_IDA = list()
  related_variable_final = NULL # correlation and MSE of final results
  epi_fisher_testresult_final = NULL # epi and fisher test results
  related_variable_pvalue_opt_final = NULL #pvalue opitimize cov matrix
  epi_test_result_mediator_detection_ByGraphStructure = NULL # find the mediator by using the graph structure, instead of IDA&JIDA.
  list_fisher_result = list()
  ######################################################################################
  #Parameter setting
  no_sim = 0 
  total_sim = 10
  try_res = 0
  print("This is next loop") # if there is something wrong with previous loop, go next
  
  ######################################################################################
  #Get the senario info.
  snps_num = MyData[x,1]
  genes_num = MyData[x,2]
  sample_num = MyData[x,3] # Sample_Total
  sample_num_Getx = MyData[x,4] # sample num consistent in Gtex,# Real gene expression data.
  randomDAG_prob = MyData[x,5] # Sparse degree of the graph
  randomDAG_lB = MyData[x,6]  # the min weight in the initial random DAG
  randomDAG_uB = abs(randomDAG_lB)*1.2 # the max weight in the initial random DAG updated by fyn 24/4/2021
  TopK_genes = round(genes_num*0.1)  # set the core gene num for fisher's exact test 
  
  s = snps_num # number of snps
  g = genes_num # number of genes 
  p = 1 # number of phenotype
  node_num = s+g+p # the node num in DAG network
  sample_num = sample_num 
  sample_num2 = sample_num # updated by fyn 22/1/2020, change the sample strategye. generate sample_num sample, and then choose sample_num_Getx of them to learn the gene-gene network, which is consistent with the real situation.
  MAFmin = 0.05
  
  
  repeat{
    
    fold_path = paste0(result_path,paste0("SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,
                                          "_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"/"))
    root_fold_path = fold_path
    
    
    if (file.exists(fold_path))
    {
      fold_path = paste0(fold_path,paste0(no_sim+1,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx"
                                          ,sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"/"))
      dir.create(fold_path)
    }
    else
    {
      dir.create(fold_path)
      fold_path2 = paste0(fold_path,paste0(no_sim+1,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",
                                           sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"/"))
      dir.create(fold_path2)
      fold_path = fold_path2
    }
    

    try_res = try({
      
      #############################################################################
      # First Step: setting the parameter in rmvDAG function in R
      # reference  https://www.rdocumentation.org/packages/pcalg/versions/2.2-0/topics/rmvDAG
      # node in this network: snps num:s /gene num: g/ phenotype num:p
      # setting weight=0 if the relation don't exist in the real world
      #############################################################################
      
      n1 = s+1
      n2 = s+g
      n3 = s+g+p
      all_serial = seq(1,n2)  #n2 is the number of gene and snps
      
      ######################################################################################
      # Get and save snp_index and gene_index
      gene_index = sample(n2,g,replace = FALSE) #sample(x, size, replace = FALSE, prob = NULL)
      snp_index = all_serial[!all_serial %in% gene_index] # after select gene index randomly, the remain part of the index is belong to snps.
      
      gene_index = sort(gene_index) #unsure about the influence of the order, maybe gene_index need to be sort to ensure there were no infulence caused by the order
      snp_index = sort(snp_index)
      
      save(gene_index, file=paste0(fold_path,"gene_index.Rdata"))
      save(snp_index, file=paste0(fold_path,"snp_index.Rdata"))
      
      
      ######################################################################################
      # Simulation the overall network
      adj = Simulate_Network(s,g,p,node_num,all_serial,gene_index,snp_index) # amat[i,j] = 0; amat[j,i]=1 represent that there is an edge i->j in amat
      adj_filepath = paste0(fold_path,paste0("SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB),
                            "_True_Graph_adj_",Sys.time(),".Rdata")
      save(adj, file=adj_filepath)


      ######################################################################################
      # Simulation the noise, then generate samples
      eMat = Simulate_Noise(sample_num2,s,g)
      rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
      rDAG = igraph.to.graphNEL(rDAG)
      samples_data2 <- rmvDAG(sample_num2, rDAG, errMat = eMat)

      
      ######################################################################################
      # Separate the sample data.
      all_serial = seq(1,sample_num2)
      sample_num_predict_imputedgene = sample_num-sample_num_Getx

      sample_index_row  = sample(nrow(samples_data2),sample_num_predict_imputedgene)

      save(sample_index_row, file=paste0(fold_path,"sample_index_row.Rdata"))

      samples_data = samples_data2[sample_index_row,] # use for prediction
      cross_validation_index = all_serial[!all_serial %in% sample_index_row]
      cross_validation_data = samples_data2[cross_validation_index,]
      
      save(samples_data2,file = paste0(fold_path,"samples_data2.Rdata"))
      save(sample_num_predict_imputedgene,file = paste0(fold_path,"sample_num_predict_imputedgene.Rdata"))
      save(cross_validation_index,file = paste0(fold_path,"cross_validation_index.Rdata"))
      index_gene2gene_samples_data = sample(nrow(samples_data2),sample_num_Getx)
      gene2gene_samples_data = samples_data2[index_gene2gene_samples_data,]
      samples_data_bootstrap = gene2gene_samples_data[,gene_index]
      save(index_gene2gene_samples_data,file = paste0(fold_path,"index_gene2gene_samples_data.Rdata"))
      
      
      ######################################################################################
      # GeneNetwork Inference
      result_list = bootstrap_pc_network3(samples_data_bootstrap,times_bootstrap,cut_off)
      graph_est = result_list[[1]]
      save(graph_est,file = paste0(fold_path,"PC_estGraph_Bootstrap.Rdata"))
      
      
      ######################################################################################
      # Get true graph, est graph and random graph, for further comparing.
      adj2 = adj[gene_index,gene_index] # we need gene network only
      rDAG1 = graph.adjacency(adj2, mode="directed", weighted=TRUE)
      rDAG1 = igraph.to.graphNEL(rDAG1)
      save(rDAG1,file = paste0(fold_path,"rDAG1.Rdata"))
      
      g1 = rDAG1 # true graph
      g2 = graph_est
      g3 <- randomDAG(g, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)
      
      
      ######################################################################################
      # Impute Gene expression data first.
      samples_for_impute = samples_data
      samples_data_iptg = Imupte_GeneExpressionData(samples_data,cross_validation_data)
      save(samples_data_iptg, file=paste0(fold_path,paste0("SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB),
                                          "_Inputed_Gene_data ",Sys.time(),".Rdata"))
      
      
      ######################################################################################
      # Caculate the Correlation between raw gene and imputed gene.
      dt_genedata = Caculate_Cor_ImputeGeneAndRawGene(gene_index,samples_data,samples_data_iptg,sample_num)
            correlation_filename = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,
                                    "_Correlation_","TopK_genes",TopK_genes,"_",Sys.time(),".csv")
      write.csv(dt_genedata, file = correlation_filename,row.names=FALSE)
      

      
      ######################################################################################
      # Gene Phenotype Network Inference
      pcSimple.fit = GeneToPhenotype_Network(samples_data_iptg)
      file_pcsimple = paste0( fold_path,"PCSimple ", "SNP", snps_num,"_Gene",genes_num,"_Sample",sample_num,"_","TopK_genes",TopK_genes,"_","randomDAG_prob_", randomDAG_prob, "randomDAG_lB_", randomDAG_lB, 
                              "pcSimple.fit_", Sys.time(), ".Rdata")	
      save(pcSimple.fit, file=file_pcsimple)
      
      
      ######################################################################################
      # Merge Gene-Gene and Gene-Phenotype Network
      bgadd2 = MergeGetx_UKBB(pcSimple.fit,graph_est)
      
      
      ######################################################################################
      # Get True Graph and Est Graph
      GenePhenotype_graph_index = union(gene_index,c(node_num))
      ture_rDAG = adj[GenePhenotype_graph_index,GenePhenotype_graph_index] #ture_rDAG the same as std_matrix
      tureMatrix = ture_rDAG
      
      estimateMatrix = bgadd2  
      estimateMatrix = t(estimateMatrix)
   
      tmp_num = g+1 # obtain some gene and phenotype data from estimateMatrix
      estimateMatrix = cbind(estimateMatrix[,2:tmp_num],estimateMatrix[,1:1])
      estimateMatrix = rbind2(estimateMatrix[2:tmp_num,],estimateMatrix[1:1,])
      
      colnames(estimateMatrix)[tmp_num] = tmp_num
      rownames(estimateMatrix)[tmp_num] = tmp_num
      
      pc.fit = as(estimateMatrix, "graphNEL")
    
      estimateGraph = graph.adjacency(estimateMatrix, mode="directed", weighted=TRUE)
      estimateGraph = igraph.to.graphNEL(estimateGraph)
      
      tureGraph = graph.adjacency(tureMatrix, mode="directed", weighted=TRUE)
      tureGraph = igraph.to.graphNEL(tureGraph)
      as(tureGraph,"matrix")

      TrueGraph_filename = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_",
                                  "TrueGarph_",Sys.time(),".Rdata")
      EstimateGraph_filename = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_",
                                      "EstimateGarph_Bootstrap_",Sys.time(),".Rdata")
      
      save(tureGraph, file= TrueGraph_filename)
      save(estimateGraph, file= EstimateGraph_filename)
      
      
      ######################################################################################
      # Caculate the mediator Info according to graph structure.
      GenePhenotype_graph_index = union(gene_index,c(node_num))
      adj_gene_outcome = adj[GenePhenotype_graph_index,GenePhenotype_graph_index] 
      mediator_result_true = Check_Mediator2(adj_gene_outcome)
      ifmediator_gene_true = mediator_result_true[[1]] # for the differ=0, we need to further analyze if it's mediator or not.the value is the mediator num
      gene_type_true = mediator_result_true[[2]] # gene_type, 2 for both direct&indirect, 0 for only mediator, 1 for only direct. to outcome, 3 for non
      
      mediator_result_est = Check_Mediator2(estimateMatrix)
      ifmediator_gene_est = mediator_result_est[[1]] # for the differ=0, we need to further analyze if it's mediator or not.the value is the mediator num
      gene_type_est = mediator_result_est[[2]] # gene_type, 2 for both direct&indirect, 0 for only mediator, 1 for only direct. to outcome, 3 for non
      
      mediator_result_true_filename = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_",
                                             "mediator_result_true_",Sys.time(),".Rdata")
      mediator_result_est_filename = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_",
                                            "mediator_result_est_",Sys.time(),".Rdata")

      save(mediator_result_true, file= mediator_result_true_filename)
      save(mediator_result_est, file= mediator_result_est_filename)
      
      
      ######################################################################################
      # Caculate the mediator Info according to graph structure.
      trueMatrix_adj = Delete_Weight_For_adj(tureMatrix)
      tab <- table(estimateMatrix,trueMatrix_adj)[2:1,2:1]
      whole_graph_epi_result = epi.tests(tab, conf.level = 0.95)
      wholeGraph_epi_result_filename = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_","wholeGraph_epi_result_",Sys.time(),".Rdata")
      save(whole_graph_epi_result, file= wholeGraph_epi_result_filename)
      
      
      ######################################################################################
      # Get the univariate results.
      univariate_test_Matrix = Univariate_Test(gene_index,samples_data,samples_data_iptg)
   
 
      ######################################################################################
      # Try Different Strategy for p-value adjustment, and then summarize the overall performance.
      n2= s+g+p
      submatrix_gene = samples_data[,gene_index] # gene data
      submatrix_phenotype = samples_data[,n2] # phenotype data
      dat1 = cbind(submatrix_gene,submatrix_phenotype) #raw gene expression data
      
      n2= s+g+p
      submatrix_gene2 = samples_data_iptg[,gene_index] # gene data
      submatrix_phenotype2 = samples_data_iptg[,n2] # phenotype data
      dat2 = cbind(submatrix_gene2,submatrix_phenotype2) #impute gene expression data
      
      
      #try to convert graphNEL to igraph in order to obtain the edgelist.
      nodes_tureGraph = nodes(tureGraph)
      phenotype_loc = length(nodes_tureGraph)  # phenotype is the last node in the Graph
      
      dt_IDA = data.frame()
      max_IDAlength1 = 0

      cov_imputedExpr = cov(dat2)
      cov_rawExpr = cov(dat1)
      
      pValue_rho = GetPvalueMatrix_ForAdjustment(result_list,pcSimple.fit,g,p)
      cov_imputedExpr2 = glassoFast(cov_imputedExpr, rho=pValue_rho)$w
      
      # after iteratively adjust by p-value from gene-gene network.
      cov_imputedExpr5 = Pvalue_adjustment_Iterative(iterateNum = 20,thresh = 0.0001,cov_imputedExpr)
      
      
      
      #method1: this is the original cov_imputedExpr without any modification
      dt_IDA1 = Get_dtIDA2(tureGraph,cov_imputedExpr,cov_imputedExpr,result_list,pcSimple.fit,genes_num,pc.fit)
      Vx2 = grep("true_jointIDA", colnames(dt_IDA1))+1
      colnames(dt_IDA1)[Vx2] = "est_jointIDA"
      
      dt_IDA1_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", 
                            "dt_IDA1_",Sys.time(), ".csv")
      write.csv(dt_IDA, file = dt_IDA1_path , row.names=TRUE)
      
      #method2: use pvalue to optimize the cov_imputedExpr, get cov_imputedExpr2
      dt_IDA2 = Get_dtIDA2(tureGraph,cov_imputedExpr2,cov_imputedExpr2,result_list,pcSimple.fit,genes_num,pc.fit)
      Vx2 = grep("true_jointIDA", colnames(dt_IDA2))+1
      colnames(dt_IDA2)[Vx2] = "est_jointIDA"
      
      dt_IDA2_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", 
                            "dt_IDA2_",Sys.time(), ".csv")
      write.csv(dt_IDA2, file = dt_IDA2_path , row.names=TRUE)
      
      
      #method5: use pvalue regularize cov_imputedExpr iteratively. 
      # important: Please be reminded that this method only iterate cov_imputedExpr regularization 
      # by pvalue, without iterate network for pc&pc_select. plese see how to get cov_imputedExpr5 for detail.
      dt_IDA5 = Get_dtIDA2(tureGraph,cov_imputedExpr5,cov_imputedExpr5,result_list,pcSimple.fit,genes_num,pc.fit)
      Vx2 = grep("true_jointIDA", colnames(dt_IDA5))+1
      colnames(dt_IDA5)[Vx2] = "est_jointIDA"
      
      dt_IDA5_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", 
                            "dt_IDA5_",Sys.time(), ".csv")
      write.csv(dt_IDA5, file = dt_IDA5_path , row.names=TRUE)
      
      
      ######################################################################################
      # Detail Analysis for dt_IDA (Initial data)
      
      #add some rank inforamtion
      dt_IDA = dt_IDA1
      ncol_IDA = ncol(dt_IDA)
      order_true_IDA = order(abs(dt_IDA$true_IDA),decreasing = TRUE)
      order_estimate_IDA = order(abs(dt_IDA$Estimate_IDA),decreasing = TRUE)
      
      
      dt_IDA[,ncol_IDA+1] = order_true_IDA
      dt_IDA[,ncol_IDA+2] = order_estimate_IDA
      
      colnames(dt_IDA)[ncol_IDA+1] = "Rank_ABS_true_IDA"
      colnames(dt_IDA)[ncol_IDA+2] = "Rank_ABS_estimate_IDA"
      
      top_gene_trueIDA_list = dt_IDA$Rank_ABS_true_IDA[1:TopK_genes] # obtain topKgenes in trueIDA list
      top_gene_estimateIDA_list = dt_IDA$Rank_ABS_estimate_IDA[1:TopK_genes] #obtain topK genes in estimateIDA list
      intersection = intersect(top_gene_trueIDA_list,top_gene_estimateIDA_list) # obtain the intersection of these two list, and the series stand the gene serial number
      
      
      contingency_table_11 = length(intersection)
      contingency_table_12 = TopK_genes-contingency_table_11
      contingency_table_21 = TopK_genes-contingency_table_11
      contingency_table_22 = genes_num-contingency_table_21-contingency_table_11-contingency_table_12
      
      fisher_test_result = fisher.test(rbind(c(contingency_table_11,contingency_table_12),c(contingency_table_21,contingency_table_22)), alternative="greater")$p.value
      
      dt_IDA[1,ncol(dt_IDA)+1] = fisher_test_result
      colnames(dt_IDA)[ncol(dt_IDA)] = "fisher_test_result"
      
     
      # save the result of univariate test into dt_IDA
      dt_IDA[,ncol(dt_IDA)+1] = univariate_test_Matrix[,1]
      colnames(dt_IDA)[ncol(dt_IDA)] = "Gene_Num"
      dt_IDA[,ncol(dt_IDA)+1] = univariate_test_Matrix[,2]
      colnames(dt_IDA)[ncol(dt_IDA)] = "Gene_rowdata_regression_Cof"
      dt_IDA[,ncol(dt_IDA)+1] = univariate_test_Matrix[,3]
      colnames(dt_IDA)[ncol(dt_IDA)] = "Gene_imputedata_regression_Cof"
      dt_IDA[,ncol(dt_IDA)+1] = univariate_test_Matrix[,4]
      colnames(dt_IDA)[ncol(dt_IDA)] = "Gene_rowdata_regression_pval"
      dt_IDA[,ncol(dt_IDA)+1] = univariate_test_Matrix[,5]
      colnames(dt_IDA)[ncol(dt_IDA)] = "Gene_imputedata_regression_pval"
      
      # order the regression coefficient in different cases
      order_Gene_rowdata_regression_Cof = order(abs(dt_IDA$Gene_rowdata_regression_Cof),decreasing = TRUE)
      dt_IDA[,ncol(dt_IDA)+1] = order_Gene_rowdata_regression_Cof
      colnames(dt_IDA)[ncol(dt_IDA)] = "Rank_Gene_rowdata_regression_Cof"
      
      
      order_Gene_imputedata_regression_Cof = order(abs(dt_IDA$Gene_imputedata_regression_Cof),decreasing = TRUE)
      dt_IDA[,ncol(dt_IDA)+1] = order_Gene_imputedata_regression_Cof
      colnames(dt_IDA)[ncol(dt_IDA)] = "Rank_Gene_imputedata_regression_Cof"
      
      
      top_gene_rowdata_regression_Cof = order_Gene_rowdata_regression_Cof[1:TopK_genes]
      top_gene_imputedata_regression_Cof = order_Gene_imputedata_regression_Cof[1:TopK_genes]
      
      
      # proportion of TopK gene estimateIDA in TopK gene trueIDA
      TopK_tureIDA_estIDA = length(intersect(top_gene_trueIDA_list,top_gene_estimateIDA_list))/length(top_gene_trueIDA_list)
      TopK_tureIDA_rowdata = length(intersect(top_gene_trueIDA_list,top_gene_rowdata_regression_Cof))/length(top_gene_trueIDA_list)
      TopK_tureIDA_imputedata = length(intersect(top_gene_trueIDA_list,top_gene_imputedata_regression_Cof))/length(top_gene_trueIDA_list)
      
      
      dt_IDA[1,ncol(dt_IDA)+1] = TopK_tureIDA_estIDA
      colnames(dt_IDA)[ncol(dt_IDA)] = "Proportion_TopK_tureIDA_estIDA"
      
      dt_IDA[1,ncol(dt_IDA)+1] = TopK_tureIDA_rowdata
      colnames(dt_IDA)[ncol(dt_IDA)] = "Proportion_TopK_tureIDA_rowdata"
      
      dt_IDA[1,ncol(dt_IDA)+1] = TopK_tureIDA_imputedata
      colnames(dt_IDA)[ncol(dt_IDA)] = "Proportion_TopK_tureIDA_imputedata"

    
    }) #end of try loop
    
    
    if (class(try_res)!="try-error") 
    {
      no_sim = no_sim + 1; 
      list_dt_IDA[[no_sim]] = dt_IDA
      
      ######################################################################################
      # Using temporarily
      dt_IDA = GetMoreDetail_OverallPerformance(dt_IDA)
      dt_IDA_final = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_","TopK_genes",TopK_genes,"_","Sim", no_sim+1, "_", 
                            "New_dt_IDA_",Sys.time(), ".csv")
      write.csv(dt_IDA, file = dt_IDA_final, row.names=FALSE)
     
      
      ######################################################################################
      # Get Fisher Test Results
      fisher_IDA_impExpr = Fisher_IDA_impExpr(dt_IDA)
      fisher_IDA = Fisher_IDA(dt_IDA)
      fisher_jointIDA_impExpr = Fisher_jointIDA_impExpr(dt_IDA)
      fisher_jointIDA = Fisher_jointIDA(dt_IDA)
      fisher_result = rbind(c("fisher_IDA_impExpr","",""),fisher_IDA_impExpr,c("fisher_IDA","",""),fisher_IDA,c("fisher_jointIDA_impExpr","",""),fisher_jointIDA_impExpr,c("fisher_jointIDA","",""),fisher_jointIDA)
      
      fisher_result_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_","TopK_genes",TopK_genes,"_", "Sim", no_sim+1, "_", "Fisher_result_",Sys.time(), ".csv")
      write.csv(fisher_result, file = fisher_result_path, row.names=FALSE)
      list_fisher_result[[no_sim]] = fisher_result
      
      
      ######################################################################################
      # Find the mediator in the gene-outcome graph by using graph structure.
      list_mediator_rval_combine = Performance_mediatorDetect_FromGraphStructure(ifmediator_gene_true,ifmediator_gene_est,genes_num)
      mediator_rval_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_","TopK_genes",TopK_genes,"_", "Sim", no_sim+5, "_", 
                                  "epi.test_mediator_FromGraphStructure_",Sys.time(), ".csv")
      write.csv(list_mediator_rval_combine, file = mediator_rval_path, row.names=FALSE)
      
      mediator_rval_combine = list_mediator_rval_combine[["mediator_rval_combine"]]
      mediator_rval_combine_colname = list_mediator_rval_combine[["mediator_rval_combine_colname"]]
      
      
      ######################################################################################
      # Get the PC_Simple performance Results
      PCsimple_mat = PC_Simple_Performance(adj,gene_index,n3,pcSimple.fi)
      PCsimple_mat_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_","TopK_genes",TopK_genes,"_","randomDAG_prob_", randomDAG_prob, 
                                 "randomDAG_lB_", randomDAG_prob, "randomDAG_uB", randomDAG_uB, "_", "Sim", no_sim+1, "_", 
                                 "PCsimple_mat_",Sys.time(), ".Rdata")
      save(PCsimple_mat, file = PCsimple_mat_path)
      
      
      ######################################################################################
      # Get the PC_Simple performance Results
      Pheno_vec = adj[gene_index,n3]
      ImpExpr_mat = Impute_Performance(dt_IDA,MyData,Pheno_vec)
      ImpExpr_mat_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_","TopK_genes",TopK_genes,"_","randomDAG_prob_", randomDAG_prob, 
                                "randomDAG_lB_", randomDAG_prob, # the min weight in the initial random DAG
                                "randomDAG_uB", randomDAG_uB, "_", "Sim", no_sim+1, "_", "ImpExpr_mat_",Sys.time(), ".Rdata")
      save(ImpExpr_mat, file = ImpExpr_mat_path)
      
      
      ######################################################################################
      # Get the summary info
      list_epi_fisher_testresult = Get_epi_fisher_Summary(PCsimple_mat,ImpExpr_mat)
      epi_fisher_testresult_colname = list_epi_fisher_testresult[["epi_fisher_testresult_colname"]]
      epi_fisher_testresult = list_epi_fisher_testresult[["epi_fisher_testresult"]]
      
      list_related_variable = Get_related_variable(dt_IDA,adj,g1,g2,g3,gene_index)
      related_variable = list_related_variable[["related_variable"]]
      related_variable_colname = list_related_variable[["related_variable_colname"]]
      
      list_related_variable_pvalue_opt = Get_related_variable_pvalue_opt(dt_IDA,dt_IDA2,dt_IDA5)
      related_variable_pvalue_opt = list_related_variable_pvalue_opt[["related_variable_pvalue_opt"]]
      related_variable_pvalue_opt_colname = list_related_variable_pvalue_opt[["related_variable_pvalue_opt_colname"]]
      
      ######################################################################################
      #correlation between the actual gene expression and outcome (should not be too high) 
      cor(dat1[,1:g], dat1[,(g+p)]) 
      cor_gene_outcome = cor(dat1[,1:g], dat1[,(g+p)]) 
      cor_gene_outcome_path = paste0(fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_","cor_gene_outcome_",Sys.time(), ".csv")
      write.csv(cor_gene_outcome, file = cor_gene_outcome_path, row.names=FALSE)
      
      
      ######################################################################################
      #Re-Organize the results. 
      if(no_sim==1)
      {
        related_variable_final = rbind(related_variable_colname,related_variable)
        epi_fisher_testresult_final = rbind(epi_fisher_testresult_colname,epi_fisher_testresult)
        related_variable_pvalue_opt_final = rbind(related_variable_pvalue_opt_colname,related_variable_pvalue_opt)
        epi_test_result_mediator_detection_ByGraphStructure = rbind(mediator_rval_combine_colname,mediator_rval_combine)
      }
      # else if(no_sim >1)
      else if(no_sim >1)
      {
        related_variable_final = rbind(related_variable_final,related_variable)
        epi_fisher_testresult_final = rbind(epi_fisher_testresult_final,epi_fisher_testresult)
        related_variable_pvalue_opt_final = rbind(related_variable_pvalue_opt_final,related_variable_pvalue_opt)
        epi_test_result_mediator_detection_ByGraphStructure = rbind(epi_test_result_mediator_detection_ByGraphStructure,mediator_rval_combine)
      }
    } 
    # debug used
    else
    {
      unlink(fold_path, recursive = TRUE) # if there some error,just delete it.
    }
    if (no_sim >= total_sim) {break}
  } #end of repeat loop
  
  ######################################################################################
  # Get the average and median for all the scenarios.
  related_variable_average = NULL
  related_variable_median = NULL
  
  ######################################################################################
  # Get the average and median for related_variable_final
  for(j in 1:ncol(related_variable_final))
  {
    average = mean(as.numeric(related_variable_final[2:nrow(related_variable_final),j]), na.rm = TRUE)
    median = median(as.numeric(related_variable_final[2:nrow(related_variable_final),j]), na.rm = TRUE)
    if(j==1) 
    { 
      related_variable_average = average 
      related_variable_median = median
    }
    else
    { 
      related_variable_average = cbind(related_variable_average,average)
      related_variable_median = cbind(related_variable_median,median)
    }
  }
  related_variable_final = rbind(related_variable_final,related_variable_average,related_variable_median)
  
  colnames(related_variable_final) <- NULL
  rownames(related_variable_final) <- seq(0,nrow(related_variable_final)-1)
  rownames(related_variable_final)[nrow(related_variable_final)-1] = "Average"
  rownames(related_variable_final)[nrow(related_variable_final)] = "Median"
  related_variable_final_path = paste0(root_fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", "related_variable_final_",Sys.time(), ".csv")
  write.csv(related_variable_final, file = related_variable_final_path , row.names=TRUE)
  

  ######################################################################################
  # Get the average and median for related_variable_pvalue_opt
  related_variable_pvalue_opt_average = NULL
  related_variable_pvalue_opt_median = NULL
  for(j in 1:ncol(related_variable_pvalue_opt_final))
  {
    average = mean(as.numeric(related_variable_pvalue_opt_final[2:nrow(related_variable_pvalue_opt_final),j]), na.rm = TRUE)
    median = median(as.numeric(related_variable_pvalue_opt_final[2:nrow(related_variable_pvalue_opt_final),j]), na.rm = TRUE)
    if(j==1) 
    { 
      related_variable_pvalue_opt_average = average 
      related_variable_pvalue_opt_median = median
    }
    else
    { 
      related_variable_pvalue_opt_average = cbind(related_variable_pvalue_opt_average,average)
      related_variable_pvalue_opt_median = cbind(related_variable_pvalue_opt_median,median)
    }
  }
  related_variable_pvalue_opt_final = rbind(related_variable_pvalue_opt_final,related_variable_pvalue_opt_average,related_variable_pvalue_opt_median)
  
  colnames(related_variable_pvalue_opt_final) <- NULL
  rownames(related_variable_pvalue_opt_final) <- seq(0,nrow(related_variable_pvalue_opt_final)-1)
  rownames(related_variable_pvalue_opt_final)[nrow(related_variable_pvalue_opt_final)-1] = "Average"
  rownames(related_variable_pvalue_opt_final)[nrow(related_variable_pvalue_opt_final)] = "Median"
  related_variable_pvalue_opt_final_path = paste0(root_fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", "related_variable_pvalue_opt_final_",Sys.time(), ".csv")
  write.csv(related_variable_pvalue_opt_final, file = related_variable_pvalue_opt_final_path , row.names=TRUE)
  
  ######################################################################################
  # caculate average and median of correlation variables 
  epi_fisher_testresult_average = NULL
  epi_fisher_testresult_median = NULL
  for(j in 1:ncol(epi_fisher_testresult_final))
  {
    average = mean(as.numeric(epi_fisher_testresult_final[2:nrow(epi_fisher_testresult_final),j]), na.rm = TRUE)
    median = median(as.numeric(epi_fisher_testresult_final[2:nrow(epi_fisher_testresult_final),j]),na.rm = FALSE)
    
    if(j==1) 
    {
      epi_fisher_testresult_average = average 
      epi_fisher_testresult_median = median
    }
    else
    {
      epi_fisher_testresult_average = cbind(epi_fisher_testresult_average,average)
      epi_fisher_testresult_median = cbind(epi_fisher_testresult_median,median)
    }
  }
  epi_fisher_testresult_final = rbind(epi_fisher_testresult_final,epi_fisher_testresult_average,epi_fisher_testresult_median)
  
  colnames(epi_fisher_testresult_final) <- NULL
  rownames(epi_fisher_testresult_final) <- seq(0,nrow(epi_fisher_testresult_final)-1)
  rownames(epi_fisher_testresult_final)[nrow(epi_fisher_testresult_final)-1] = "Average"
  rownames(epi_fisher_testresult_final)[nrow(epi_fisher_testresult_final)] = "Median"
  epi_fisher_testresult_final_path = paste0(root_fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", "epi_fisher_testresult_final_",Sys.time(), ".csv")
  write.csv(epi_fisher_testresult_final, file = epi_fisher_testresult_final_path , row.names=TRUE)
  
  
  ######################################################################################
  #added by fyn 29/3/2021. caculate the average of epi.test for mediator dection from graphstructure.
  epi_test_result_mediator_detection_FromGraphstructure_average = NULL
  epi_test_result_mediator_detection_FromGraphstructure_median = NULL
  for(j in 1:ncol(epi_test_result_mediator_detection_ByGraphStructure))
  {
    average = mean(as.numeric(epi_test_result_mediator_detection_ByGraphStructure[2:nrow(epi_test_result_mediator_detection_ByGraphStructure),j]), na.rm = TRUE)
    median = median(as.numeric(epi_test_result_mediator_detection_ByGraphStructure[2:nrow(epi_test_result_mediator_detection_ByGraphStructure),j]),na.rm = FALSE)
    
    if(j==1) 
    {
      epi_test_result_mediator_detection_FromGraphstructure_average = average 
      epi_test_result_mediator_detection_FromGraphstructure_median  = median
    }
    else
    {
      epi_test_result_mediator_detection_FromGraphstructure_average = cbind(epi_test_result_mediator_detection_FromGraphstructure_average,average)
      epi_test_result_mediator_detection_FromGraphstructure_median = cbind(epi_test_result_mediator_detection_FromGraphstructure_median,median)
    }
  }
  epi_test_result_mediator_detection_ByGraphStructure = rbind(epi_test_result_mediator_detection_ByGraphStructure,epi_test_result_mediator_detection_FromGraphstructure_average,epi_test_result_mediator_detection_FromGraphstructure_median)
  colnames(epi_test_result_mediator_detection_ByGraphStructure) <- NULL
  rownames(epi_test_result_mediator_detection_ByGraphStructure) <- seq(0,nrow(epi_test_result_mediator_detection_ByGraphStructure)-1)
  rownames(epi_test_result_mediator_detection_ByGraphStructure)[nrow(epi_test_result_mediator_detection_ByGraphStructure)-1] = "Average"
  rownames(epi_test_result_mediator_detection_ByGraphStructure)[nrow(epi_test_result_mediator_detection_ByGraphStructure)] = "Median"
  mediator_detection_ByGraphStructure_path = paste0(root_fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", "epi.test_results_mediator_detection_ByGraphStructure_",Sys.time(), ".csv")
  write.csv(epi_test_result_mediator_detection_ByGraphStructure, file = mediator_detection_ByGraphStructure_path , row.names=TRUE)
  ###################################End#################################################
  
  
  # caculate average
  fisher_result_average = Fisher_Results_Average(list_fisher_result)
  colnames(fisher_result_average) <-c("topK_percent","OR","pval_IDA")
  
  fisher_result_average_path = paste0(root_fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", "fisher_result_average_",Sys.time(), ".csv")
  write.csv(fisher_result_average, file = fisher_result_average_path , row.names=FALSE)
  
  # caculate median
  fisher_result_median = Fisher_Results_Median(list_fisher_result)
  colnames(fisher_result_median) <-c("topK_percent","OR","pval_IDA")
  
  fisher_result_median_path = paste0(root_fold_path,"SNP",snps_num,"_Gene",genes_num,"_Sample",sample_num,"_sample_num_Getx",sample_num_Getx,"_randomDAG_prob",randomDAG_prob,"_randomDAG_lB",randomDAG_lB,"_", "fisher_result_median_",Sys.time(), ".csv")
  write.csv(fisher_result_median, file = fisher_result_median_path , row.names=FALSE)
  
}

proc.time() - dt1



