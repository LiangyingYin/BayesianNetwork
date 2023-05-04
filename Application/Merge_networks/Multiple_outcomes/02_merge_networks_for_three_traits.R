
#*************************************************************
# Written by: yinly, August 19,2020
# Merge graphs with 2 clinical traits as outcome
#*************************************************************
library(pcalg)
library(graph)
library(dplyr)

Get_Merged_BN_Graph <- function(whole_graph_name,whole_graph_selected_genes_name,UKBB_name,predicted_trait_name){
    # load BN graph with the "cause" trait as the outcome for the BN, e.g., BMI for BMI on DM
    # whole_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Exploratory_Analysis/Merged_Graphs_Relax_Rand/"
    whole_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Merged_Graphs_Update/"
    #whole_graph_name = "Whole_Blood_BMI_Graph.Rdata"
    load(paste0(whole_graph_prefix,whole_graph_name))
    whole_graph_matrix = bgadd # After first merge, the object name has been revised as bgadd
    whole_graph_matrix = t(whole_graph_matrix)# to make bgadd consistent with the graph derived from UKBB, [i,j]=1 indicate the direction from i to j 
    
    #load selected genes file for whole graph of "cause" trait
    # whole_graph_selected_genes_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_All_Selected_Genes/"
    whole_graph_selected_genes_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Merged_Genes_Results_Update/"
    #whole_graph_selected_genes_name = "Whole_Blood_BMI_Selected_Genes.csv" #selected genes for the "cause" trait, BMI in our BMI->CAD example
    whole_graph_selected_genes = fread(paste0(whole_graph_selected_genes_prefix,whole_graph_selected_genes_name))

    # load pcSimple.fit with imputed gene expression and genetically-predicted "cause" trait as predictors and "effect" trait as outcome
    UKBB_prefix = "/exeh_3/yinly/BayesianNetwork/03_UKBB_Network/02_Multiple_Traits/Order3_Results_CIS_Trans/"
    #UKBB_name = "UKBB_BMI_CAD_Whole_Blood.Rdata"
    load(paste0(UKBB_prefix,UKBB_name))  
    genes_names = as.vector(names(pcSimple.fit$G))
    assoc = as.numeric(as.vector(pcSimple.fit$G))
    zMin = as.vector(pcSimple.fit$zMin)
    UKBB_result = as.data.frame(cbind(genes_names,assoc,zMin))
    # Extract the causal relationship between genetically predicted traits and studied outcome, e.g., BMI and CAD
    UKBB_Traits_assoc = UKBB_result[nrow(UKBB_result),]
    
    # Merge BN graph for the "cause" trait with pcSimple.fit for all traits
    Merged_result = left_join(whole_graph_selected_genes,UKBB_result,by="genes_names")
    #predicted_trait_name = "ENSG_BMI"
    predicted_traits_outcome = data.frame(predicted_trait_name,"NA","NA","NA","NA",UKBB_Traits_assoc[2],UKBB_Traits_assoc[3])
    names_Merged_result <- c("genes_names","assoc.x","zMin.x","assoc.y","zMin.y","assoc.z","zMin.z")
    colnames(Merged_result) <- names_Merged_result
    #Merged_result <- droplevels(Merged_result)
    names(predicted_traits_outcome) = names_Merged_result
    Merged_result_final = rbind(Merged_result,predicted_traits_outcome,fill=TRUE)
    #Merged_result_final <- droplevels(Merged_result)


    Merged_result_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Merged_Genes_Results_Update/"
    Merged_result_name = gsub("Rdata","csv",UKBB_name)
    fwrite(Merged_result_final,paste0(Merged_result_prefix,Merged_result_name),sep="\t")  

    genes_to_outcome = as.vector(Merged_result$assoc.z)
    # Revised by yinly, June 18, 2021
    # Description: if directly combine with UKBB_Traits_assoc[[2]], we get 2. The revised version solve this problem.
    assoc <- as.numeric(UKBB_Traits_assoc[[2]]==1)
    predictors_to_outcome = as.vector(c(genes_to_outcome,assoc))
    predictors_to_outcome = as.numeric(predictors_to_outcome)
    # Find the index of genes with NA, then replace them with 0
    NA_index = which(is.na(predictors_to_outcome))
    predictors_to_outcome[NA_index] = 0
    # Revised by yinly: adjacency matrix derived from whole graph with the "cause" trait as the outcome was combined first
    final_mat = cbind(whole_graph_matrix,predictors_to_outcome) ##[i,j] means an edge from i to j (here it is gene -> disease, but NOT the other way round
    final_mat <- rbind(final_mat,rep(0,(ncol(whole_graph_matrix)+1))) ##need to check if the diagonals are all zero 
    #final_mat = as.matrix(final_mat)
    class(final_mat) = "numeric"

    #*******************************************************
    # Merge two graph by add background knowledge
    #*******************************************************
    #the gInput takes [i,j] to mean an edge from j to i, which is opposite to the convention in graph package
    # Commented by:yinly, Date: August 6, 2020. If we have convert the adjacency matrix to graphNEL object, there is no need to transpose the adjacency matrix
    mat_input = t(final_mat) 
   
    t1 = proc.time()

    bgadd = addBgKnowledge(gInput = mat_input,  
                        verbose = TRUE, 
                        checkInput = FALSE)					   
  
    proc.time()-t1

    Merged_graph_prefix = "/exeh_3/yinly/BayesianNetwork/04_Merge_GTEx_UKBB/Multiple_Traits/Merged_Graphs_Update/"
    Merged_graph_name = gsub("UKBB_","",UKBB_name)
    save(bgadd,file=paste0(Merged_graph_prefix,Merged_graph_name))
}

 info_mat = as.matrix(fread("Multi_Traits_Graph_Merge_Dict_Update_1April.csv",header=TRUE))
 index <-c(9,11)
 for(j in 1:length(index)){
 #for(i in 8:11){
    i <- index[j]
    whole_graph_name = info_mat[i,1]
    whole_graph_selected_genes_name = info_mat[i,2]
    UKBB_name = info_mat[i,3]
    predicted_trait_name = info_mat[i,4]
    Get_Merged_BN_Graph(whole_graph_name,whole_graph_selected_genes_name,UKBB_name,predicted_trait_name)
    cat("The graph merge for",gsub(".Rdata","",UKBB_name),"is completed!\n")
}
#whole_graph_name = "BMI_Diabetes_Whole_Blood.Rdata"
#whole_graph_selected_genes_name = "UKBB_BMI_Diabetes_Whole_Blood.csv"
#UKBB_name = "UKBB_BMI_Diabetes_Response_CAD_Whole_Blood.Rdata"
#predicted_trait_name = "ENSG_Diabetes"
#Get_Merged_BN_Graph(whole_graph_name,whole_graph_selected_genes_name,UKBB_name,predicted_trait_name)
#cat("The graph merge for",gsub(".Rdata","",UKBB_name),"is completed!\n")
