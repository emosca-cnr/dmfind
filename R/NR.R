#' Network Resampling
#' @param G igraph object
#' @param rankedVector ranked named vector with vertex scores; only `names(rankedVector)` in `V(G)$name` will be considered
#' @param k number of permutations
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using `SerialParam()`. See ?bplapply.
#' @return list of values
#' @export
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom igraph induced.subgraph get.adjacency V
#' @importFrom stats setNames

NR <- function(G = NULL, rankedVector = NULL, k = 99, BPPARAM = NULL) {
    
    if (is.null(BPPARAM)) {
      BPPARAM <- SerialParam()
    }
    
    cat("BPPARAM\n")
    print(BPPARAM)
    
    stopifnot(all(names(rankedVector) %in% V(G)$name))
    
    #real values
    rankedVectorNorm <- rankedVector / max(rankedVector)
    
    G <- induced.subgraph(G, V(G)$name[match(names(rankedVectorNorm), V(G)$name)])
    omegaVect <- omega(G = G, u = rankedVectorNorm)
    
    #permutations of rankedVectorNorm: named vector of indices
    idxPerm <- lapply(1:k, function(x)  setNames(sample(1:length(rankedVectorNorm)), names(rankedVectorNorm)))
    
    Gi <- induced.subgraph(G, V(G)$name[match(names(rankedVectorNorm), V(G)$name)])
    Ai <- as.matrix(get.adjacency(Gi))
    res <- bplapply(idxPerm, omega_perm, dS = rankedVectorNorm, Ai = Ai, BPPARAM = BPPARAM)
    
    #n-by-k matrix
    res <- t(do.call(rbind, res))
    
    #LMC p values
    out <- apply(res, 2, function(x) sign(x >= omegaVect))
    out <- rowSums(out)
    out <- (out + 1) / (k + 1)
    
    #p_adj <- p.adjust(out, method = "fdr")
    #q_val <- eFDR(real_values = omegaVect, all_values = c(omegaVect, as.numeric(res)), BPPARAM = BPPARAM)
    
    return(list(
      NRsummary = data.frame(
        id = names(out),
        rank = 1:length(rankedVectorNorm),
        rankingScore = rankedVectorNorm,
        omega = omegaVect,
        p_val = out,
        stringsAsFactors = FALSE
      ),
      omegaPerm = res
    ))
    
  }
