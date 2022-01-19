#' Calculate ccf
#' 
#' @description Calculate ccf
#' @param G igraph object
#' @param rankedVectorNames ranked vector of names
#' @return cff
#' @usage calculate_ccf(G, rankedVectorNames)
#' @examples 
#' \dontrun{calculate_ccf(G, rankedVectorNames)}
#' @import igraph
#' @export

calculate_ccf <- function(G=NULL, rankedVectorNames=NULL){

    n <- length(rankedVectorNames)
    rankedVectorNames <- 
        rankedVectorNames[1:min(length(rankedVectorNames), 1000)]
    n <- length(rankedVectorNames)

    ans <- sapply(1:n, function(i) 
        CCF(igraph::induced_subgraph(G, 
        V(G)$name[V(G)$name %in% rankedVectorNames[1:i]])))

    return(ans)

}
