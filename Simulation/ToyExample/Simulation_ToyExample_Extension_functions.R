
Simulate_Noise<-function(sample_num2,s,g)
{
  set.seed(123)
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

Imupte_GeneExpressionData<-function(samples_data,cross_validation_data)
{
  # allocate snps data and gene data randomly
  set.seed(123)
  snps_data = cross_validation_data[,snp_index]  ##the 800 subjects (GTEx)
  gene_data = cross_validation_data[,gene_index] ##the 800 subjects (GTEx)
  snps_data2 = samples_data[,snp_index] #use for prediction
  
  
  samples_data_iptg  = samples_for_impute
  
  
  for (k in gene_index)
  {
    
    gene_k = cross_validation_data[,k] # obtain the gene data in turn 
    
    cvfit =cv.glmnet(snps_data,gene_k,nfolds = cv.glmnet_nfolds)
    
    coef(cvfit, s = "lambda.min")
    
    
    #choosing the best lambda according to cross-validation by using  cv.glmnet
    pred_genes = predict(cvfit, newx = snps_data2[1:sample_num_predict_imputedgene,], s = "lambda.min")
    if (var(pred_genes)==0) {
      
      cvm_min_nonzero = min( cvfit$cvm[cvfit$nzero>0] )  #extract the best cv error at which there is at least one non-zero coeff. 
      lamb = cvfit$lambda[cvfit$cvm==cvm_min_nonzero]
      pred_genes = predict(cvfit, newx = snps_data2[1:sample_num_predict_imputedgene,], s = lamb)
    }
    
    
    samples_data_iptg[,k] = pred_genes[,1:1]  #use the 1st column data temporarily
    
    
  }
  
  return(samples_data_iptg)
  
}


bootstrap_pc_network3<-function(samples_data_bootstrap,times_bootstrap,cut_off)
{
  set.seed(123)
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
    # if(i==1) #or randomly choose
    # {
    #   #set.seed 
    #   save(new_samples, file=paste0(fold_path,"new_samples.Rdata"))
    # }
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

Caculate_Cor_ImputeGeneAndRawGene<-function(gene_index,samples_data,samples_data_iptg,sample_num)
{
  set.seed(123)
  dt_genedata = data.frame()
  row_count = 1
  
  for (g_count in gene_index) # iterate over each gene according to gene_index.
  {
    
    v1 = samples_data[,g_count]
    
    v2 = samples_data_iptg[,g_count]
    
    # results in comparsion
    MSE_value = MSE(v1,v2)/sample_num 
    cor_value = cor(v1,v2)  # in some cases, cor value is NA, because the standard deviation is zero
    cosine_value = cosine(v1,v2) 
    
    dt_genedata[row_count,1] = paste0("G",g_count)
    dt_genedata[row_count,2] = MSE_value
    dt_genedata[row_count,3] = cor_value
    dt_genedata[row_count,4] = cosine_value
    colnames(dt_genedata)[1] = "Gene_num"
    colnames(dt_genedata)[2] = "MSE_value"
    colnames(dt_genedata)[3] = "cor_value"
    colnames(dt_genedata)[4] = "cosine_value"
    row_count = row_count+1
  }
  
  return(dt_genedata) #print gene imputed data and simulated data
  
}

GeneToPhenotype_Network<-function(samples_data_iptg)
{
  
  set.seed(123)
  library(dplyr) 
  library(data.table)
  
  var_thres = 0.01
  target_outcome_name = "X20016.0.0"  #Diastolic blood pressure, automated reading
  
  pheno_dir = "/exeh_4/yaning_feng/04_Simulation/Data/array_phenos.txt"
  
  column3 = s+g+1
  column4 = s+g+p
  
  pheno_data = as.matrix(samples_data_iptg[1:sample_num_predict_imputedgene, column3:column4])
  
  
  gene_data = samples_data_iptg[1:sample_num_predict_imputedgene,gene_index] # extract the gene data according to gene_index
  
  
  dfnew_rankNorm = cbind(pheno_data,gene_data) #鐠嬪啯宕瞘ene閸滃henotype閻ㄥ嫭鏆熼幑顕嗙礉閹跺henotype閺佺増宓侀弨鎯ф躬閸撳秹�??
  colnames(dfnew_rankNorm)[1] <-"outcome"
  
  
  library(coop) 
  
  expr_corMat = pcor(dfnew_rankNorm[,-1])
  
  
  cor_with_outcome = apply(dfnew_rankNorm[,-1],  2 , function(col) {pcor(dfnew_rankNorm[,"outcome"], col)}  )
  precompute_corMat = cbind(cor_with_outcome, expr_corMat)
  
  precompute_corMat = rbind( c(1,cor_with_outcome), precompute_corMat) 
  
  
  #********************************************************
  # Apply PC-simple algorithm around the response variable on NFBC
  #**************************************************************
  library(pcalg)
  #library(kpcalg)
  
  ## Load  data
  n <- nrow (dfnew_rankNorm)
  V <- colnames(dfnew_rankNorm) # labels aka node names
  
  t1 = proc.time()
  ## estimate local causal network structure around the response variable 
  
  
  pcSimple.fit <- pcSelect_stableVer(y=dfnew_rankNorm[,1], 
                                     dm=dfnew_rankNorm[,-1] , 
                                     corMat = precompute_corMat, 
                                     alpha=0.001, 
                                     max_ord=3, 
                                     corMethod = "standard",
                                     verbose = TRUE, directed = TRUE)  
  
  proc.time()-t1
  file_path = "/exeh/exe4/yaning_feng/BayesianNetwork/04_Simulation/Results/fyntest_simulate.Rdata"
  # save(pcSimple.fit, file=file_path)
  return(pcSimple.fit)
}



MergeGetx_UKBB<-function(pcSimple.fit,graph_est)
{
  set.seed(123)
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

Univariate_Test<-function(gene_index,samples_data,samples_data_iptg)
{
  set.seed(123)
  outcome = samples_data[,dim(samples_data)[2]] # phenotype or outcome is the last column in this sample matrix
  
  
  univariate_test_Matrix = matrix(nrow = length(gene_index), ncol = 5)
  
  count = 1
  for(gene in gene_index)
  {
    print(gene)
    gene_rowdata = samples_data[,gene]
    gene_imputedata = samples_data_iptg[,gene]
    
    lm_gene_rowdata = lm(outcome~gene_rowdata)
    lm_gene_imputedata = lm(outcome~gene_imputedata)
    
    lm_gene_rowdata_pval = summary(lm_gene_rowdata)$coefficients[2,4]
    lm_gene_imputedata_pval = summary(lm_gene_imputedata)$coefficients[2,4]
    
    # univariate_test_Matrix[count,1] = gene # write gene_index into this column 
    univariate_test_Matrix[count,1] = count
    univariate_test_Matrix[count,2] = unname(coefficients(lm_gene_rowdata)[2])
    univariate_test_Matrix[count,3] = unname(coefficients(lm_gene_imputedata)[2])
    univariate_test_Matrix[count,4] = lm_gene_rowdata_pval
    univariate_test_Matrix[count,5] = lm_gene_imputedata_pval
    
    count = count + 1
  }
  return(univariate_test_Matrix)
}

GetPvalueMatrix_ForAdjustment<-function(result_list,pcSimple.fit,g,p)
{
  set.seed(123)
  pValue_gg = result_list[[2]]
  
  adj_gg = as(result_list[[1]],"matrix") #update by 9/6/2021.
  pValue_gg[which(adj_gg==0 & pValue_gg < 0.01, arr.ind = T)] = 1 #update by 9/6/2021. Here 0.01 is the alpha from PC-algorithm.
  
  pValue_gp = 2*pnorm(-abs(pcSimple.fit$zMin))
  pValue_rho = matrix(nrow= g+p, ncol=g+p)
  pValue_rho [1:g,1:g] = pValue_gg
  pValue_rho [,g+p] <- c( pValue_gp,-1)
  pValue_rho [g+p,] <- c( pValue_gp,-1)
  
  diag(pValue_rho) <- 0
  pValue_rho[pValue_rho==-1] <- 1
  return(pValue_rho)
}

Pvalue_adjustment_Iterative<-function(iterateNum = 20,thresh = 0.0001,cov_imputedExpr,pValue_rho)
{
  set.seed(123)
  iterate_error2 = array(0, iterateNum)
  cov_imputedExpri2 = cov_imputedExpr
  
  for(i in 1:iterateNum)
  {
    GLASSO = glassoFast(cov_imputedExpri2, rho=pValue_rho)
    cov_imputedExpri2 = GLASSO$w
    iterate_error2[i] = iterate_error2[i] + sum(GLASSO$wi[upper.tri(GLASSO$wi)] != 0) *
      log(nrow(dat2))/2
    if((i>1)&&(abs(iterate_error2[i]-iterate_error2[i-1])< thresh))
    {
      break
    }
    
  }
  
  cov_imputedExpr5 = cov_imputedExpri2
  return(cov_imputedExpr5)
}

Get_dtIDA2 <- function(tureGraph,cov_rowExpr,cov_imputedExpr,result_list,pcSimple.fit,genes_num,pc.fit)
{
  set.seed(123)
  #try to convert graphNEL to igraph in order to obtain the edgelist.
  nodes_tureGraph = nodes(tureGraph)
  phenotype_loc = length(nodes_tureGraph)  # phenotype is the last node in the Graph
  nodes(tureGraph) = as.character(seq(1,phenotype_loc))
  nodes_tureGraph = nodes(tureGraph)
  
  dt_IDA = data.frame()
  max_IDAlength1 = 0
  
  for(gene in 1:genes_num)
  {
    node1 = nodes_tureGraph[gene] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    
    covTrue <- trueCov(tureGraph)
    true_IDA <- ida(node1,node2, covTrue, tureGraph, method = "local", type = "pdag")
    
    l.IDA.estimate <- ida(node1,node2, cov_imputedExpr , pc.fit, method = "local", type = "pdag")							   
    
    
    print(l.IDA.estimate) # print the estimate IDA value
    
    num_local_IDA = length(l.IDA.estimate)
    if(num_local_IDA > max_IDAlength1)
    {
      max_IDAlength1 = num_local_IDA
    }
  }
  
  max_IDAlength = max_IDAlength1
  #tureGraph as standard
  
  
  nodes_true = nodes(tureGraph)
  x.pos = nodes_true[!nodes_true %in% node2]
  y.pos = node2
  
  x.pos = as.numeric(x.pos) # change the data type from character to numeric
  y.pos = as.numeric(y.pos)
  
  #jointIDA compute
  true_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,covTrue,covTrue,graphEst=tureGraph,technique="RRC")
  
  #estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov(dat2),graphEst=estimateGraph,technique="RRC")
  estimate_jointIDA_list = jointIda(x.pos=x.pos,y.pos=y.pos,cov_rowExpr,cov_imputedExpr,graphEst = pc.fit,technique="RRC")
  
  
  for(gene in 1:genes_num)
  {
    node1 = nodes_tureGraph[gene] # node1 is gene
    node2 = phenotype_loc # node2 is phenotype in this case
    
    #IDA compute
    covTrue <- trueCov(tureGraph)
    true_IDA <- ida(node1,node2, covTrue, tureGraph, method = "local", type = "pdag")
    # true_IDA = causalEffect(tureGraph,node2,node1)
    
    
    l.IDA.estimate <- ida(node1,node2, cov_imputedExpr , pc.fit, method = "local", type = "pdag")	
    
    jointIDA_position = which(x.pos==node1)
    true_jointIDA = true_jointIDA_list[jointIDA_position]
    
    
    dt_IDA[gene,1] = node1 # first node in the edge
    dt_IDA[gene,2] = node2 # second node in the edge
    
    dt_IDA[gene,3] = true_IDA
    num_local_IDA = length(l.IDA.estimate)
    
    
    for(numlIDA in 1:num_local_IDA)
    {
      dt_IDA[gene,(3+numlIDA)] = l.IDA.estimate[numlIDA]
    }
    
    #filled missing value
    tmp1=  4+max_IDAlength
    tmp2 = 3+numlIDA
    minus_1_2 = tmp1-tmp2
    if(minus_1_2>1)
    {
      for(tmp in 1:minus_1_2)
      {
        # dt_IDA[gene_num,(tmp2+tmp-1)] = 'NA'
        dt_IDA[gene,(tmp2+tmp)] = 'NA'
      }
    }
    
    
    dt_IDA[gene,(4+max_IDAlength)] = true_jointIDA
    
    num_estimate_jointIDA = ncol(estimate_jointIDA_list)
    for(num_jointIDA in 1:num_estimate_jointIDA)
    {
      dt_IDA[gene,(4+max_IDAlength+num_jointIDA)] = estimate_jointIDA_list[jointIDA_position,num_jointIDA]
    }
  }
  
  colnames(dt_IDA)[1] = "gene_loc"
  colnames(dt_IDA)[2] = "phenotype_loc"
  colnames(dt_IDA)[3] = "true_IDA"
  colnames(dt_IDA)[4] = "Estimate_IDA"
  colnames(dt_IDA)[(4+max_IDAlength)] = "true_jointIDA"
  
  
  print("This is dt_IDA")
  print(dt_IDA)
}

#pcSelect_stableVer function defination
library(coop)
library(pcalg)
pcSelect_stableVer <- function(y, dm, alpha, corMat = NA, corMethod = "standard", max_ord=3,
                               verbose = FALSE, directed = FALSE)
{
  ## Purpose: Find columns in dm, that have nonzero parcor with y given
  ## any other set of columns in dm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - y: Response Vector (length(y)=nrow(dm))
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - verbose: 0-no output, 1-small output, 2-details
  ## ----------------------------------------------------------------------
  ## Value: List
  ## - G: boolean vector with connected nodes
  ## - zMin: Minimal z values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 27.4.07
  
  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  ## cl <- match.call()
  
  zMin <- c(0,rep.int(Inf,p))
  #***************************
  # This line changed so we can use a faster function to calculate correlation matrix
  #****************************
  #C <- mcor(cbind(y,dm), method = corMethod)
  if (!is.na(corMat)) {C <- corMat} else { 
    C <- coop::pcor(cbind(y,dm) ) } 
  
  
  
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- c(FALSE,rep.int(TRUE,p))
  seq_p <- seq_len(p+1L) # = 1:(p+1)
  
  done <- FALSE
  ord <- 0
  
  while (!done && any(G) && ord<=max_ord) {
    G_order_specific <- G    ##G_order_specific stores which nodes are connected to the outcome for each order K 
    
    
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G_order_specific)  #which(G)
    remEdges <- length(ind)
    if(verbose >= 1)
      cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
    
    
    for (i in 1:remEdges) {
      if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=",i,"|iMax=",remEdges,"\n")
      y <- 1
      x <- ind[i]
      
      if (G[x]) {
        # nbrsBool <- G
        # nbrsBool[x] <- FALSE
        # nbrs <- seq_p[nbrsBool]
        
        #*****************************************************
        # NEW part: consider all nodes left over from previous order  as covariates 
        #******************************************************
        nbrsBool_orderSpecific <- G_order_specific ##G is the current "active" set of x that are connected to y FROM THE PREVIOUS ORDER of tests (eg if now you are testing order 1 dependency, we keep all genes significant at order 0 tests as possible 'neighbours')  
        nbrsBool_orderSpecific [x] <- FALSE  ## the predictor (xk) that you're now testing is of course not considered as 'neighbour' of itself
        nbrs_orderSpecific <- seq_p[nbrsBool_orderSpecific]  ## neighbors of y without itself and x
        
        length_nbrs <- length(nbrs_orderSpecific)
        
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)
          
          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat {
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            
            z <- zStat(x,y, nbrs_orderSpecific[S], C,n)
            if(abs(z) < zMin[x]) zMin[x] <- abs(z)
            
            if (verbose >= 2)
              cat(paste("x:",vNms[x-1],"y:",(ytmp <- round((p+1)/2)),"S:"),
                  c(ytmp,vNms)[nbrs_orderSpecific[S]], paste("z:",z,"\n"))
            
            if (abs(z) <= cutoff) {
              G[x] <- FALSE
              break
            }
            
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
              
            }
          } ## {repeat}
        }
      } ## end if( G )
    } ## end for(i ..)
    
    
    ord <- ord+1
  } ## end while
  
  ## return
  list(G = setNames(G[-1L], vNms),
       zMin = zMin[-1])
}## pcSelect


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


