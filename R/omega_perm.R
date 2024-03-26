#' Calculation of permuted omega function
#' @param dS scores vector; it must have the same names
#' and size of the vertices of G
#' @param idx vector of randomly resampled dS scores;
#' the names must be the same as dS's
#' @return vector of permuted omega
#' @import igraph


omega_perm <- function(idx = NULL, dS = NULL, Ai=NULL) {
    
    #subnetwork among the considered nodes
    #Gi <- induced.subgraph(G, V(G)$name[match(names(dS), V(G)$name)])
    #Ai <- as.matrix(get.adjacency(Gi))
    
    #n-1 lists: from 2:n elements of idx sorted according to random indices idx
    idx2n <- lapply(2:length(dS), function(x) sort(idx[names(idx) %in% names(dS)[1:x]]))
    
    
    #cross product of correct values
    omegaVect <- dS %*% t(dS)
    
    #from 2:n calculate omega
    omegaVect <-
        c(0, unlist(
            lapply(
                idx2n,
                calc_omega_i,
                Ai = Ai,
                dSprod = omegaVect
            )
        ))
    
    return(omegaVect)
    
}
