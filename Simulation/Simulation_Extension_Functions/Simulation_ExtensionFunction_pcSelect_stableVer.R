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

