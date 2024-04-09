#' Calculation of omega function
#' 
#' @description Calculation of omega function
#' @param G igraph object; only the vertices in names(u) will be considered
#' @param u ranked named list; names must correspond to V(G)$name
#' @return omega 
#' @importFrom igraph induced.subgraph get.adjacency V
#' @export

omega <- function(G=NULL, u=NULL) {

    #implementation 2
    # note that the vertex order is not the same as in u
    Gi <- induced.subgraph(G, match(names(u), V(G)$name)) 
    # extract local adjacency matrix
    Ai <- as.matrix(get.adjacency(Gi)) 
    
    # appropriate setting of Ai and u names
    idxNorm <- match(names(u), rownames(Ai)) 
    Ai <- Ai[idxNorm, idxNorm]

    U <- matrix(u, ncol = 1, dimnames = list(names(u)))
    omegaVect <-  U %*% t(U)
    omegaVect <- omegaVect * Ai
    omegaVect[upper.tri(omegaVect)] <- 0

    omegaVect <- cumsum(rowSums(omegaVect))

    return(omegaVect)
}
