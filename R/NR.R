#' Local Monte Carlo p values
#' @param G gene x gene undirected interaction graph
#' @param rankedVector scores vector; it must have the same names and size of the vertices of G
#' @param k number of permutations
#' @param mc.cores number of cores
#' @param nullModel null model type
#' @return list of values
#' @usage NR(G, rankedVector, k=99, mc.cores=1, nullModel="A", norm="n")
#' @example NR(G, rankedVector, k=99, mc.cores=1, nullModel="A", norm="n")
#' @export
#' @import BiocParallel


NR <- function(G, rankedVector, k=99, mc.cores=1, nullModel="A", norm="n") {

  nullModel <- match.arg(nullModel)

  #real values
  rankedVectorNorm <- rankedVector / max(rankedVector)

  omegaVect <- omega(G, rankedVectorNorm, norm = norm)

  #permutations of rankedVectorNorm: named vector of indices
  idxPerm <- lapply(1:k, function(x) array(sample(1:length(rankedVectorNorm)), 
             dimnames = list(names(rankedVectorNorm))))

  if(mc.cores > 1){
    res <- BiocParallel::bplapply(idxPerm, omega_perm, G=G, 
                                  dS=rankedVectorNorm, mc.cores = mc.cores, 
                                  nullModel=nullModel, norm=norm)
  }else{
    res <- lapply(idxPerm, omega_perm, G=G, dS=rankedVectorNorm, 
                  nullModel=nullModel, norm=norm)
  }

  #n-by-k matrix
  res <- t(do.call(rbind, res))

  #LMC p values
  out <- apply(res, 2, function(x) sign(x >= omegaVect))
  out <- rowSums(out)
  out <- (out + 1) / (k+1)

  return(list(NRsummary=data.frame(id=names(out), 
                                   rank=1:length(rankedVectorNorm), 
                                   rankingScore=rankedVectorNorm, 
                                   omega=omegaVect, p=out, 
                                   stringsAsFactors = FALSE), omegaPerm=res))

}


