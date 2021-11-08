#' Calculate ccf
#' @export

calculate_ccf <- function(G=NULL, rankedVectorNames=NULL){

    n <- length(rankedVectorNames)
    rankedVectorNames <- 
        rankedVectorNames[1:min(length(rankedVectorNames), 1000)]
    n <- length(rankedVectorNames)

    ans <- sapply(1:n, function(i) 
        CCF(igraph::induced_subgraph(g, 
        V(g)$name[V(g)$name %in% rankedVectorNames[1:i]])))

    return(ans)

}
