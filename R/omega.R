#' Calculation of omega function
#' 
#' @description Calculation of omega function
#' @param G igraph object; only the vertices in names(u) will be considered
#' @param u ranked named list; names must correspond to V(G)$name
#' @param norm normalize by number of links (l), number or vertices 
#' (v) or do not normalize (n)
#' @return omega 
#' @usage omega(G, u, norm=c("n", "l", "v"))
#' @examples 
#' \dontrun{omega(G, u, norm=c("n", "l", "v"))}
#' @import igraph
#' @export

omega <- function(G, u, norm=c("n", "l", "v")) {

    norm <- match.arg(norm, c("n", "l", "v"))

    #implementation 2
    # note that the vertex order is not the same as in u
    Gi <- igraph::induced.subgraph(G, match(names(u), V(G)$name)) 
    # extract local adjacency matrix
    Ai <- as.matrix(get.adjacency(Gi)) 
    # appropriate setting of Ai and u names
    idxNorm <- match(names(u), rownames(Ai)) 
    Ai <- Ai[idxNorm, idxNorm]

    U <- matrix(u, ncol = 1, dimnames = list(names(u)))
    omegaVect <-  U %*% t(U)
    omegaVect <- omegaVect * Ai
    omegaVect[upper.tri(omegaVect)] <- 0

    #omegaVect <- 2*cumsum(rowSums(omegaVect))
    omegaVect <- cumsum(rowSums(omegaVect))

    if(norm=="l"){
        Ai[upper.tri(Ai)] <- 0
        nE <- cumsum(rowSums(Ai))
        omegaVect <- omegaVect / nE
        omegaVect[is.nan(omegaVect)] <- 0
    }

    if(norm=="v"){
        omegaVect <- omegaVect / (1:length(omegaVect))
    }

    return(omegaVect)
}
