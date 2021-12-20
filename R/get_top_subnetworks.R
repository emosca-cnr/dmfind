#' Calculation of permutation-adjusted network smoothing index
#' 
#' @description Calculation of permutation-adjusted network smoothing index
#' @param G igraph graph object
#' @param rankedList ranked vector; names must match V(G)$name
#' @param idx vector of ranks according to which networks will be extracted
#' @return list with snG and snSummary elements
#' @examples get_top_subnetworks(G, rankedList, idx=1:length(rankedList))
#' @usage get_top_subnetworks(G, rankedList, idx=1:length(rankedList))
#' @import igraph
#' @export
#'

get_top_subnetworks <- function(G, rankedList, idx=1:length(rankedList)){

    if(is.null(idx)){
        idxCurr <- 1:length(rankedList)
    }else{
        idxCurr <- idx
    }
    gOut <- vector('list', length = length(idxCurr))
    statOut <- data.frame(idx=idxCurr, nV=idxCurr, nE=0, d=0, 
                           nComp=0, nVmax=0, nVmin=0, nVconnTot=0)
    j <- 1
    for(i in idxCurr){
        gOut[[j]] <- 
            induced_subgraph(G, 
                             V(G)$name[V(G)$name %in% names(rankedList)[1:i]])
        j <- j+1
    }

    statOut$nE <- unlist(lapply(gOut, function(x) length(E(x))))
    statOut$d <- unlist(lapply(gOut, graph.density))

    cl <- lapply(gOut, clusters)
    statOut$nComp <- unlist(lapply(cl, function(x) x$no))
    statOut$nVmax <- unlist(lapply(cl, function(x) max(x$csize)))
    statOut$nVmin <- unlist(lapply(cl, function(x) min(x$csize)))
    statOut$nVconnTot <- 
        unlist(lapply(cl, function(x) sum(x$csize[x$csize>1])))
    statOut$omega <- 
        unlist(lapply(idxCurr, function(x) 
            omega(G, rankedList[1:x] / max(rankedList[1:x]))[x]))

    return(list(snG=gOut, snSummary=statOut))
    
}
