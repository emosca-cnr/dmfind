#' Calculation of permutation-adjusted network smoothing index
#'
#' @param X0 input matrix
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, see normalize_adj_mat function
#' @param k number of permutations
#' @param eps numeric value
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @param returnPerm return permutations
#' @param mode whether to use network smoothing index (S) (default) or
#' network smoothing values (Xs)
#' @param ... additional parameters of ND
#' @description Calculation of permutation-adjusted network smoothing index
#' @usage calc_adjND(X0, W, eps, k=99, mode=c("S", "Xs"), returnPerm=FALSE, ...)
#' @examples
#' \dontrun{calc_adjND(X0, W, eps, k=99, mode=c("S", "Xs"), returnPerm=FALSE)}
#' @return \code{data.frame} with X0, Xs, S, p and Sp
#' @import BiocParallel
#' @export
#'
calc_adjND <-
  function(X0=NULL,
           W=NULL,
           eps = NULL,
           k = 99,
           mode = c("S", "Xs"),
           BPPARAM = NULL,
           returnPerm = FALSE,
           ...) {
    mode <- match.arg(mode)
    
    if (is.null(BPPARAM)) {
      BPPARAM <- SerialParam()
    }
    
    cat("BPPARAM\n")
    print(BPPARAM)
    
    if (is.null(eps)) {
      cat("eps value(s) not provided. It is important to tune it! Using 1...\n")
      eps <- matrix(1, nrow = 1, ncol = ncol(X0))
    }
    
    if (!is.matrix(eps)) {
      stop("eps must be a 1-row matrix with same number of columns of X0\n")
    }
    
    if (ncol(eps) != ncol(X0)) {
      stop("eps must have the same number of columns of X0\n")
    }
    
    allX0 <-
      c(list(X0), lapply(1:k, function(x)
        matrix(
          as.numeric(X0),
          ncol = ncol(X0),
          dimnames =
            list(sample(rownames(X0), nrow(X0)))
        )))
    
    allX0 <- BiocParallel::bplapply(allX0, function(x)
      x[match(rownames(W), rownames(x)), , drop = FALSE], BPPARAM = BPPARAM)
    
    cat("network propagation\n")
    
    allXS <- BiocParallel::bplapply(allX0, function(x)
      ND(x, W, ...)$Xs, BPPARAM = BPPARAM)
    
    
    if (mode == "S") {
      cat("calculation of Sp\n")
      #S
      all_S <-
        lapply(1:(k + 1), function(i)
          nsi(allX0[[i]], allXS[[i]], eps = eps))
      Xs <- allXS[[1]]
    } else{
      all_S <- allXS
    }
    
    ## trash all permutations
    if (!returnPerm) {
      rm(allX0, allXS)
    }
    
    #p
    estP <- calc_p(all_S)
    S <- all_S[[1]]
    
    if (mode == "S") {
      if (returnPerm) {
        out <- list(
          X0 = X0,
          Xs = Xs,
          eps = eps,
          S = S,
          p = estP,
          Sp = S * -log10(estP),
          allX0 = allX0,
          allXS = allXS
        )
      } else{
        out <- list(
          X0 = X0,
          Xs = Xs,
          eps = eps,
          S = S,
          p = estP,
          Sp = S * -log10(estP)
        )
      }
    } else{
      if (returnPerm) {
        out <- list(
          X0 = X0,
          Xs = S,
          eps = eps,
          p = estP,
          Xsp = S * -log10(estP),
          allX0 = allX0,
          allXS = allXS
        )
      } else{
        out <- list(
          X0 = X0,
          Xs = S,
          eps = eps,
          p = estP,
          Xsp = S * -log10(estP)
        )
      }
    }
    
    return(out)
  }
