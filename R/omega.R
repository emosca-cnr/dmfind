#' Calculation of omega function
#' 
#' @description Calculation of omega function
#' @param G igraph object; only the vertices in names(u) will be considered
#' @param u ranked named list; names must correspond to V(G)$name
#' @param norm normalize by number of links (l), number or vertices (v) or do not normalize (n)
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
    idx.norm <- match(names(u), rownames(Ai)) 
    Ai <- Ai[idx.norm, idx.norm]

    U <- matrix(u, ncol = 1, dimnames = list(names(u)))
    omega_vect <-  U %*% t(U)
    omega_vect <- omega_vect * Ai
    omega_vect[upper.tri(omega_vect)] <- 0

    #omega_vect <- 2*cumsum(rowSums(omega_vect))
    omega_vect <- cumsum(rowSums(omega_vect))

    if(norm=="l"){
        Ai[upper.tri(Ai)] <- 0
        nE <- cumsum(rowSums(Ai))
        omega_vect <- omega_vect / nE
        omega_vect[is.nan(omega_vect)] <- 0
    }

    if(norm=="v"){
        omega_vect <- omega_vect / (1:length(omega_vect))
    }

    return(omega_vect)
}
