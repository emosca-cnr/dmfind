#' Calculate the connected component fraction
#' 
#' @description Calculate the connected component fraction
#' @param G igraph object
#' @param rankedVectorNames ranked vector of names
#' @param min_size min size of a connected component
#' @return cff
#' @importFrom igraph V induced_subgraph
#' @export

calculate_ccf <- function(G=NULL, rankedVectorNames=NULL, min_size=2){

    n <- length(rankedVectorNames)
    rankedVectorNames <- 
        rankedVectorNames[1:min(length(rankedVectorNames), 1000)]
    n <- length(rankedVectorNames)

    ans <- sapply(1:n, function(i) CCF(induced_subgraph(G, V(G)$name[V(G)$name %in% rankedVectorNames[1:i]]), min_size=min_size))

    return(ans)

}
