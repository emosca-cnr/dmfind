#' Calculation of permuted omega function
#' @param G gene x gene undirected interaction graph
#' @param dS scores vector; it must have the same names and size of the vertices of G
#' @param idx vector of randomly resampled dS scores; the names must be the same as dS's
#' @param null_model null model type to be used
#' @return vector of permuted omega
#' @usage omega_perm(idx, G, dS, null_model="A", norm="n")
#' @examples omega_perm(idx, G, dS, null_model="A", norm="n")
#' @import igraph


omega_perm <- function(idx, G, dS, null_model="A", norm="n") {


    #subnetwork among the considered nodes
    Gi <- igraph::induced.subgraph(G, match(names(dS), V(G)$name))
    Ai <- as.matrix(get.adjacency(Gi))

    #n-1 lists: from 2:n elements of idx sorted according to random indices idx
    idx2n <- lapply(2:length(dS), function(x) 
        sort(idx[names(idx) %in% names(dS)[1:x]]))

    #  if(null_model=="A"){

    #cross product of correct values
    omegaVect <- dS %*% t(dS)

    #from 2:n calculate omega
    omegaVect <- c(0, unlist(lapply(idx2n, calc_omega_i, 
                                     Ai=Ai, dSprod=omegaVect, norm=norm)))

    #  }
    
    # if(null_model=="Af"){
    #
    #   Aperm <- Ai[idx, idx]
    #   colnames(Aperm) <- rownames(Aperm) <- rownames(Ai)
    #   Gperm <- graph.adjacency(Aperm, mode = "undirected")
    #   #cross product of correct values
    #   omegaVect <- omega(Gperm, u = dS)
    #
    # }
    #
    #
    # if(null_model=="S"){
    #
    #   #matrix dS*dS' * Ai with random dS
    #   omegaVect <- (dS[idx] %*% t(dS[idx])) * Ai
    #   #progressive omega
    #   omegaVect <- sapply(1:length(dS), function(i) sum(omegaVect[1:i, 1:i]))
    #
    # }
    # if(null_model=="SA"){
    #
    #   #matrix dS*dS' * Ai with random dS
    #   omegaVect <- (dS[idx] %*% t(dS[idx])) * Ai[idx, idx]
    #   #progressive omega
    #   omegaVect <- sapply(1:length(dS), function(i) sum(omegaVect[1:i, 1:i]))
    #
    # }
    
    return(omegaVect)

}
