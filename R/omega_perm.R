#' Calculation of permuted omega function
#' @param G gene x gene undirected interaction graph
#' @param dS scores vector; it must have the same names and size of the vertices of G
#' @param idx vector of randomly resampled dS scores; the names must be the same as dS's
#' @param null_model null model type to be used
#' @import igraph


omega_perm <- function(idx, G, dS, null_model="A", norm=FALSE) {


  #subnetwork among the considered nodes
  Gi <- igraph::induced.subgraph(G, match(names(dS), V(G)$name))
  Ai <- as.matrix(get.adjacency(Gi))

  #n-1 lists: from 2:n elements of idx sorted according to random indices idx
  idx_2_n <- lapply(2:length(dS), function(x) sort(idx[names(idx) %in% names(dS)[1:x]]))

  #  if(null_model=="A"){

  #cross product of correct values
  omega_vect <- dS %*% t(dS)

  #from 2:n calculate omega
  omega_vect <- c(0, unlist(lapply(idx_2_n, calc_omega_i, Ai=Ai, dSprod=omega_vect, norm=norm)))

  #  }

  # if(null_model=="Af"){
  #
  #   Aperm <- Ai[idx, idx]
  #   colnames(Aperm) <- rownames(Aperm) <- rownames(Ai)
  #   Gperm <- graph.adjacency(Aperm, mode = "undirected")
  #
  #   #cross product of correct values
  #   omega_vect <- omega(Gperm, u = dS)
  #
  # }
  #
  #
  # if(null_model=="S"){
  #
  #   #matrix dS*dS' * Ai with random dS
  #   omega_vect <- (dS[idx] %*% t(dS[idx])) * Ai
  #   #progressive omega
  #   omega_vect <- sapply(1:length(dS), function(i) sum(omega_vect[1:i, 1:i]))
  #
  # }
  #
  # if(null_model=="SA"){
  #
  #   #matrix dS*dS' * Ai with random dS
  #   omega_vect <- (dS[idx] %*% t(dS[idx])) * Ai[idx, idx]
  #   #progressive omega
  #   omega_vect <- sapply(1:length(dS), function(i) sum(omega_vect[1:i, 1:i]))
  #
  # }


  return(omega_vect)

}
