######################################################################################
# parameter setting for this ToyExample: 
# SNP700_Gene70_Sample100000_sample_num_Getx8000_randomDAG_prob0.05_randomDAG_lB-0.3

######################################################################################
# Load the source code and function code
# /opt/R-3.6.1/bin/R
ROOT_FOLDER = "/exeh_4/yaning_feng/04_Simulation/Github_code/BackupV1.1/"
FOLDER_PATH = paste0(ROOT_FOLDER,"ToyExample/") # Please change this path to your local path for loading the related functions successfully

# Please change the random DAG if needed for testing, and the gene_index, snp_index should be changed accordingly
load(paste0(FOLDER_PATH,"/randomDAGs/randomDAG_sim.Rdata")) 
load(paste0(FOLDER_PATH,"/randomDAGs/gene_index.Rdata")) 
load(paste0(FOLDER_PATH,"/randomDAGs/snp_index.Rdata")) 


source(paste0(FOLDER_PATH,"Simulation_ToyExample_Extension_functions.R"))
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
result_path = FOLDER_PATH

min_weight_gene_to_pheno = 0
max_weight_gene_to_pheno = 3
max_weight_snp_to_gene = 0.5 #in real situation, relations between snps and genes is weak  
prop_snp_related_gene = 0.01  #set the proportion of snps related to gene 
prop_snp_predict_gene = 1

cv.glmnet_nfolds = 5  # for cross-validation
eMat_g_tmp_rnorm_sd = 2  #gene expression data error variance
prop_direct_causal_genes = 1/5

#some parameters about bootstrap
times_bootstrap = 1 
cut_off = 0.6


######################################################################################
# Please noted that since we load the randomDAG for saveing time, please do NOT change these settings in this toy example, in order to be conssitent with the parameter for generating the DAG. 
snps_num = 700
genes_num = 70
sample_num = 100000
sample_num_Getx = 800
randomDAG_prob = 0.05
randomDAG_lB = -0.3
randomDAG_uB = abs(randomDAG_lB)*1.2 

s = snps_num # number of snps
g = genes_num # number of genes 
p = 1 # number of phenotype
node_num = s+g+p # the node num in DAG network
sample_num = sample_num 
sample_num2 = sample_num 
MAFmin = 0.05

n1 = s+1
n2 = s+g
n3 = s+g+p

eMat = Simulate_Noise(sample_num2,s,g)
rDAG = graph.adjacency(adj, mode="directed", weighted=TRUE)
rDAG = igraph.to.graphNEL(rDAG)
samples_data2 <- rmvDAG(sample_num2, rDAG, errMat = eMat)


######################################################################################
# Separate the sample data.
all_serial = seq(1,sample_num2)
sample_num_predict_imputedgene = sample_num-sample_num_Getx

sample_index_row  = sample(nrow(samples_data2),sample_num_predict_imputedgene)
samples_data = samples_data2[sample_index_row,] # use for prediction
cross_validation_index = all_serial[!all_serial %in% sample_index_row]
cross_validation_data = samples_data2[cross_validation_index,]


index_gene2gene_samples_data = sample(nrow(samples_data2),sample_num_Getx)
gene2gene_samples_data = samples_data2[index_gene2gene_samples_data,]
samples_data_bootstrap = gene2gene_samples_data[,gene_index]


######################################################################################
# GeneNetwork Inference
result_list = bootstrap_pc_network3(samples_data_bootstrap,times_bootstrap,cut_off)
graph_est = result_list[[1]]


######################################################################################
# Get true graph, est graph and random graph, for further comparing.
adj2 = adj[gene_index,gene_index] # we need gene network only
rDAG1 = graph.adjacency(adj2, mode="directed", weighted=TRUE)
rDAG1 = igraph.to.graphNEL(rDAG1)


g1 = rDAG1 # true graph
g2 = graph_est
g3 <- randomDAG(g, prob = randomDAG_prob, lB = randomDAG_lB, uB=randomDAG_uB)


######################################################################################
# Impute Gene expression data first.
samples_for_impute = samples_data
samples_data_iptg = Imupte_GeneExpressionData(samples_data,cross_validation_data)


######################################################################################
# Caculate the Correlation between raw gene and imputed gene.
dt_genedata = Caculate_Cor_ImputeGeneAndRawGene(gene_index,samples_data,samples_data_iptg,sample_num)


######################################################################################
# Gene Phenotype Network Inference
pcSimple.fit = GeneToPhenotype_Network(samples_data_iptg)


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
cov_imputedExpr5 = Pvalue_adjustment_Iterative(iterateNum = 20,thresh = 0.0001,cov_imputedExpr,pValue_rho)


fold_path = FOLDER_PATH
#method1: this is the original cov_imputedExpr without any modification
dt_IDA1 = Get_dtIDA2(tureGraph,cov_imputedExpr,cov_imputedExpr,result_list,pcSimple.fit,genes_num,pc.fit)
Vx2 = grep("true_jointIDA", colnames(dt_IDA1))+1
colnames(dt_IDA1)[Vx2] = "est_jointIDA"


#method2: use pvalue to optimize the cov_imputedExpr, get cov_imputedExpr2
dt_IDA2 = Get_dtIDA2(tureGraph,cov_imputedExpr2,cov_imputedExpr2,result_list,pcSimple.fit,genes_num,pc.fit)
Vx2 = grep("true_jointIDA", colnames(dt_IDA2))+1
colnames(dt_IDA2)[Vx2] = "est_jointIDA"



#method5: use pvalue regularize cov_imputedExpr iteratively. 
# important: Please be reminded that this method only iterate cov_imputedExpr regularization 
# by pvalue, without iterate network for pc&pc_select. plese see how to get cov_imputedExpr5 for detail.
dt_IDA5 = Get_dtIDA2(tureGraph,cov_imputedExpr5,cov_imputedExpr5,result_list,pcSimple.fit,genes_num,pc.fit)
Vx2 = grep("true_jointIDA", colnames(dt_IDA5))+1
colnames(dt_IDA5)[Vx2] = "est_jointIDA"


Cor_TureIDA_EstIDA_Pearson = cor(dt_IDA1$true_IDA,dt_IDA1$Estimate_IDA,method = c("pearson"))
Cor_TureIDA_EstIDA_Pearson_pvalue	= cor(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA,method = c("pearson"))
Cor_TureIDA_EstIDA_Pearson_pvalue_iterate	= cor(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA,method = c("pearson"))
Cor_TureJointIDA_EstJointIDA_Pearson = cor(dt_IDA1$true_jointIDA,dt_IDA1$est_jointIDA,method = c("pearson"))
Cor_TureJointIDA_EstJointIDA_Pearson_pvalue	= cor(dt_IDA2$true_jointIDA,dt_IDA2$est_jointIDA,method = c("pearson"))
Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate	= cor(dt_IDA5$true_jointIDA,dt_IDA5$est_jointIDA,method = c("pearson"))
RMSE_TrueIDA_EstIDA	= RMSE(dt_IDA1$true_IDA,dt_IDA1$Estimate_IDA)
RMSE_TrueIDA_EstIDA_pvalue	= RMSE(dt_IDA2$true_IDA,dt_IDA2$Estimate_IDA)
RMSE_TrueIDA_EstIDA_pvalue_iterate	= RMSE(dt_IDA5$true_IDA,dt_IDA5$Estimate_IDA)
RMSE_TrueJointIDA_EstJointIDA	= RMSE(dt_IDA1$true_jointIDA,dt_IDA1$est_jointIDA)
RMSE_TrueJointIDA_EstJointIDA_pvalue = RMSE(dt_IDA2$true_jointIDA,dt_IDA2$est_jointIDA)	
RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate = RMSE(dt_IDA5$true_jointIDA,dt_IDA5$est_jointIDA)	

comp_res = c(Cor_TureIDA_EstIDA_Pearson,
             Cor_TureIDA_EstIDA_Pearson_pvalue,	
             Cor_TureIDA_EstIDA_Pearson_pvalue_iterate,	
             Cor_TureJointIDA_EstJointIDA_Pearson,
             Cor_TureJointIDA_EstJointIDA_Pearson_pvalue,	
             Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate,	
             RMSE_TrueIDA_EstIDA,	
             RMSE_TrueIDA_EstIDA_pvalue,	
             RMSE_TrueIDA_EstIDA_pvalue_iterate,
             RMSE_TrueJointIDA_EstJointIDA,	
             RMSE_TrueJointIDA_EstJointIDA_pvalue,	
             RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate)

col_names = c("Cor_TureIDA_EstIDA_Pearson",
              "Cor_TureIDA_EstIDA_Pearson_pvalue",	
              "Cor_TureIDA_EstIDA_Pearson_pvalue_iterate",	
              "Cor_TureJointIDA_EstJointIDA_Pearson",
              "Cor_TureJointIDA_EstJointIDA_Pearson_pvalue",	
              "Cor_TureJointIDA_EstJointIDA_Pearson_pvalue_iterate",	
              "RMSE_TrueIDA_EstIDA",	
              "RMSE_TrueIDA_EstIDA_pvalue",	
              "RMSE_TrueIDA_EstIDA_pvalue_iterate",
              "RMSE_TrueJointIDA_EstJointIDA",	
              "RMSE_TrueJointIDA_EstJointIDA_pvalue",	
              "RMSE_TrueJointIDA_EstJointIDA_pvalue_iterate")

dt_res = rbind(col_names,comp_res)
dt_res = as.data.table(dt_res)
fwrite(dt_res,file = paste0(FOLDER_PATH,"simulation_ToyExample_res.csv"))

