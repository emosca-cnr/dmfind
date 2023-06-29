#' Local Monte Carlo p values
#' @param G gene x gene undirected interaction graph
#' @param rankedVector scores vector; it must have the same
#' names and size of the vertices of G
#' @param k number of permutations
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @param nullModel null model type
#' @param norm type of normalization
#' @return list of values
#' @usage NR(G, rankedVector, k=99, mc.cores=1, nullModel="A", norm="n")
#' @examples
#' \dontrun{
#' NR(G, rankedVector, k=99, mc.cores=1, nullModel="A", norm="n") }
#' @export
#' @import BiocParallel
#' @import igraph


NR <-
  function(G = NULL,
           rankedVector = NULL,
           k = 99,
           nullModel = "A",
           norm = "n",
           BPPARAM = NULL) {
    nullModel <- match.arg(nullModel)
    
    if (is.null(BPPARAM)) {
      BPPARAM <- SerialParam()
    }
    
    cat("BPPARAM\n")
    print(BPPARAM)
    
    
    #real values
    rankedVectorNorm <- rankedVector / max(rankedVector)
    
    omegaVect <- omega(G, rankedVectorNorm, norm = norm)
    
    #permutations of rankedVectorNorm: named vector of indices
    idxPerm <-
      lapply(1:k, function(x)
        array(sample(1:length(rankedVectorNorm)),
              dimnames = list(names(rankedVectorNorm))))
    
    res <-
      BiocParallel::bplapply(
        idxPerm,
        omega_perm,
        G = G,
        dS = rankedVectorNorm,
        BPPARAM = BPPARAM,
        norm = norm
      )
    
    #n-by-k matrix
    res <- t(do.call(rbind, res))
    
    #LMC p values
    out <- apply(res, 2, function(x)
      sign(x >= omegaVect))
    out <- rowSums(out)
    out <- (out + 1) / (k + 1)
    
    return(list(
      NRsummary = data.frame(
        id = names(out),
        rank = 1:length(rankedVectorNorm),
        rankingScore = rankedVectorNorm,
        omega = omegaVect,
        p = out,
        stringsAsFactors = FALSE
      ),
      omegaPerm = res
    ))
    
  }
