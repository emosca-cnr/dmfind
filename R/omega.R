#' Calculation of omega function
#' @param G gene x gene undirected interaction graph
#' @param u ranked named list; names have to be those in V(G)$name
#' @import igraph


omega <- function(G, u) {

  #implementation 2
  Gi <- igraph::induced.subgraph(G, match(names(u), V(G)$name)) # note that the vertex order is not the same as in u
  Ai <- as.matrix(get.adjacency(Gi)) # extract local adjacency matrix
  idx.norm <- match(names(u), rownames(Ai)) # appropriate setting of Ai and u names
  Ai <- Ai[idx.norm, idx.norm]

  U <- matrix(u, ncol = 1, dimnames = list(names(u)))
  omega_vect <-  U %*% t(U)
  omega_vect <- omega_vect * Ai
  omega_vect[upper.tri(omega_vect)] <- 0
  omega_vect <- 2*cumsum(rowSums(omega_vect))

  return(omega_vect)
}
