#' Calculation of permutation-adjusted network smoothing index
#'
#' @param X0 input matrix
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, see normalize_adj_mat function
#' @param k number of permutations
#' @param eps matrix of eps values; same columns of X0; see eval_eps().
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(); See BiocParallel::bplapply()
#' @param bin_type see NPATools::perm_vertices()
#' @param method see NPATools::perm_vertices()
#' @param cut_par see NPATools::perm_vertices()
#' @param vertex_sets see NPATools::perm_vertices()
#' network smoothing values (Xs)
#' @param ... additional arguments passed to NPATools::ND()
#' @return list of data frames 
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom  NPATools perm_vertices perm_X0 ND calc_p eFDR
#' @importFrom  igraph degree graph_from_adjacency_matrix
#' @export
#'
calc_adjND <-
  function(X0 = NULL,
           W = NULL,
           eps = NULL,
           k = 99,
           BPPARAM = NULL,
           bin_type = "number",
           method = "simple",
           cut_par = 20,
           vertex_sets = NULL,
           ...) {

    if (is.null(BPPARAM)) {
      BPPARAM <- SerialParam()
    }
    
    cat("BPPARAM\n")
    print(BPPARAM)
    
    if (is.null(eps)) {
      stop("Missing eps value(s). It is important to tune it! See eval_eps().\n")
    }
    stopifnot(is.matrix(eps), ncol(eps) == ncol(X0))
    
    vert_deg <- degree(graph_from_adjacency_matrix(sign(W), mode = "undirected"))
    
    #vertex set permutations
    vert_perms <- perm_vertices(vert_deg = vert_deg, k=k, method = method, cut_par = cut_par, bin_type = bin_type, vertex_sets = vertex_sets)
    
    cat("ND over X0 permutations\n")
    allX0 <- perm_X0(X0 = X0, perms = vert_perms)
    allX0 <- lapply(allX0, function(x) x[match(rownames(W), rownames(x)), , drop = FALSE])
    
    allXSX0 <- bplapply(allX0, function(x) ND(X0 = x, W = W, ...), BPPARAM = BPPARAM)
    
    cat("ND over W permutations\n")
    allW <- perm_X0(X0 = W, perms = vert_perms)
    allW <- lapply(allW, function(x) x[match(rownames(W), rownames(x)), match(colnames(W), colnames(x))])
    
    allXSW <-  bplapply(allW, function(x) ND(X0 = X0, W = x, ...), BPPARAM = BPPARAM)
    
    
    #save real Xs
    Xs <- allXSX0[[1]] #allXSW[[1]] is the same
    
    cat("calculation of S\n")
    
    #S perm X0
    all_SX0 <- lapply(1:(k + 1), function(i) nsi(X0 = allX0[[i]], Xs = allXSX0[[i]], eps = eps))
    
    #S perm W
    all_SW <- lapply(1:(k + 1), function(i) nsi(X0 = X0, Xs = allXSW[[i]], eps = eps))
    
    ## trash all permutations
    rm(allX0, allXSX0, allXSW)
  
    S <- all_SX0[[1]] #true S,  same as all_SW[[1]]
    
    cat("calculation of p-values\n")
    out <- list(
      X0 = X0,
      Xs = Xs,
      eps = eps,
      S = S,
      pX0 = calc_p(all_SX0),
      pW = calc_p(all_SW)
    )
    
    out$SpX0 <- out$S * -log10(out$pX0)
    out$SpW <- out$S * -log10(out$pW)
    
    ### probability is the are under the hyperbole of the two p-values, that is:
    # p = p1p2
    # p' = p - p * ln(p)
    # this is the same as the Fisher combination based on x-squared distribution, that is:
    # poolr::fisher()
    
    out$p <- out$pX0 * out$pW
    out$p <- out$p - (out$p * log(out$p))
      
    out$Sp <- out$S * -log10(out$p)
      
    # cat("calculation of FDR\n")
    # e_fdr <- X0
    # for (i in 1:ncol(e_fdr)) {
    #   all_values <- unlist(lapply(all_SX0, function(S_perm) S_perm[, i])) #the ith-column of every permutation
    #   e_fdr[, i] <- eFDR(real_values = S[, i], all_values = all_values)
    # }
    # out$FDRX0 <- e_fdr
    # 
    # e_fdr <- X0
    # for (i in 1:ncol(e_fdr)) {
    #   all_values <- unlist(lapply(all_SW, function(S_perm) S_perm[, i])) #the ith-column of every permutation
    #   e_fdr[, i] <- eFDR(real_values = S[, i], all_values = all_values)
    # }
    # out$FDRW <- e_fdr

    return(out)
  }
