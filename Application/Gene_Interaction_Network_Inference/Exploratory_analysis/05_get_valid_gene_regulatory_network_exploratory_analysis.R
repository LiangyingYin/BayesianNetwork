
library(pcalg)
library(data.table)
library(graph)

#**************************************************************
# Written by yinly, Oct 29,2020
# Description: re-oriented undirected edges based on 
# edge information derived from CAM
#**************************************************************

#**********************************************************************************************
# Goal: Get undirected edges from pc_obj
# param@pc_obj: object derived from PC alogrithm
# param@u(return): matrix define the undirected edges
#**********************************************************************************************
get_undirected_edges <- function (object, labels = NULL) 
{
    g <- getGraph(object)
    if (is.null(labels)) 
        labels <- nodes(g)
    isDir <- isDirected(g)
    wm <- wgtMatrix(g)
    if (isDir && !is(object, "pcAlgo")) 
        stop("implementation only for 'pcAlgo', currently ...")
    wmU <- wm + t(wm)
    wmD <- t(wm - t(wm))
    u <- which(wmU == 2 & upper.tri(wmU), arr.ind = TRUE)
    return(u)
}
#**********************************************************************************************
# Goal: Get undirected edges from pc_obj
# param@pc_obj: object derived from PC alogrithm
# param@u(return): proportion of the undirected edges over all identified edges
#**********************************************************************************************
get_reorientation_proportion <-function(object, labels = NULL){
    g <- getGraph(object)
    if (is.null(labels)) 
        labels <- nodes(g)
    isDir <- isDirected(g)
    wm <- wgtMatrix(g)
    if (isDir && !is(object, "pcAlgo")) 
        stop("implementation only for 'pcAlgo', currently ...")
    wmU <- wm + t(wm)
    wmD <- t(wm - t(wm))
    u <- which(wmU == 2 & upper.tri(wmU), arr.ind = TRUE)
    cat("The number of undirected edges is",nrow(u),"\n")
    if (isDir) {
        d <- which(wmD == 1, arr.ind = TRUE)
        d <- d[order(d[, 1]), ]
        cat("The number of directed edges is",nrow(d),"\n")
        prop = nrow(u)/(nrow(u) + nrow(d))
        return(prop)
    }
    else {
        d <- matrix(0, 0, 0)
        cat("Proportion calculation is not applicable for this object","\n")
        return(NULL)
    }
}
#**********************************************************************************************
# Goal: Get directed edges from pc_obj
# param@pc_obj: object derived from PC alogrithm
# param@d(return): matrix define the directed edges
#**********************************************************************************************
get_directed_edges <- function (object, labels = NULL) 
{
    g <- getGraph(object)
    if (is.null(labels)) 
        labels <- nodes(g)
    isDir <- isDirected(g)
    wm <- wgtMatrix(g)
    if (isDir && !is(object, "pcAlgo")) 
        stop("implementation only for 'pcAlgo', currently ...")
    wmU <- wm + t(wm)
    wmD <- t(wm - t(wm))
    d <- which(wmD == 1, arr.ind = TRUE)
    d <- d[order(d[,1]),]
    return(d)
}

GTEx_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/GTEx_Network_Update/"
CAM_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/CAM_Network_Update/"
file_Dict = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/InvalidGraphs_Update.csv"
Reoriented_GTEx_prefix = "/exeh_3/yinly/BayesianNetwork/02_GTex_Network/01_Single_Tissue/Exploratory_Analysis/Reoriented_GTEx_Network_Update/"
Dict = as.matrix(fread(file_Dict,header=TRUE))
for(i in 1:nrow(Dict)){
    GTEx_name = Dict[i,2]
    load(paste0(GTEx_prefix,GTEx_name))##object name: pc_obj
    trait_name = gsub("_Causal_Network_PCSelect.Rdata","",GTEx_name)
    CAM_name = Dict[i,6]
    load(paste0(CAM_prefix,CAM_name))
    cat("Graph merge for the",i,"trait",gsub("_Causal_Network_PCSelect.Rdata","",GTEx_name),"has started!","\n")
    pc_graph = pc_obj@graph
    adj_mat = as(pc_graph, "matrix")
    #*********************************************************************************
    # Revised by:yinly, Oct 29,2020
    # Description: reorient undirected edges as directed edge based on edge information
    # derived from CAM algorithm
    #*********************************************************************************
    cam_adj = as.matrix(cam_obj$Adj) 
    proportion_undirected = get_reorientation_proportion(pc_obj)
    cat("The proportion of undirected edges is",proportion_undirected,"\n")
    pMax = pc_obj@pMax
    if(proportion_undirected>0){
       undirected_egdes = get_undirected_edges(pc_obj)
       # re-orient undirected edges as directed edges
       # Define a variable to record missing edges
       missing_edges = 0
       for(i in 1:nrow(undirected_egdes)){
           row_ind = undirected_egdes[i,1]
           col_ind = undirected_egdes[i,2]
           # Re-orient undirected based on edge informations derived from cam_obj
           adj_mat[row_ind,col_ind] = cam_adj[row_ind,col_ind]
           adj_mat[col_ind,row_ind] = cam_adj[col_ind,row_ind]
           if((cam_adj[row_ind,col_ind]==0) & (cam_adj[col_ind,row_ind]==0)){
               missing_edges = missing_edges + 1
           }
        }
        cat("There is",missing_edges,"edges missing for",trait_name,"after edge re-orientation based on CAM algorithm! \n")
        # Revised by yinly, June 13, 2021
        # Description: get ambigous v-structures(coded as 2 in adj_mat) and re-orient them.
        ambigous_edges <- which(adj_mat==2,arr.ind=TRUE)
        for(j in 1:nrow(ambigous_edges)){
            row_ind_ambig <- ambigous_edges[j,1]
            col_ind_ambig <- ambigous_edges[j,2]
            adj_mat[row_ind_ambig,col_ind_ambig] <- cam_adj[row_ind_ambig,col_ind_ambig]
        }
        # Revised by yinly, June 13, 2021
        # Description, get directed edges and redirect them based on CAM objects as well to see if it can solve cycles.
        directed_egdes <- get_directed_edges(pc_obj)
        for(k in 1:nrow(directed_egdes)){
            row_ind_d <- directed_egdes[k,1]
            col_ind_d <- directed_egdes[k,2]
            adj_mat[row_ind_d,col_ind_d] <- cam_adj[row_ind_d,col_ind_d]
        }

        graph.fit = as(adj_mat, "graphNEL")
        pc_obj@graph = graph.fit
        save(pc_obj,file= paste0(Reoriented_GTEx_prefix,GTEx_name))
        #adj_temp = as(pc_obj@graph,"matrix")
        #amat = t(adj_temp)
        #cat("This is the",i,"graph",trait_name,"\n")
        #isValidGraph(amat,verbose=TRUE)
    }
}

