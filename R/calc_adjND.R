#' Calculation of permutation-adjusted network smoothing index
#' @param X0 input matrix
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, see normalize_adj_mat function
#' @param k number of permutations
#' @param eps numeric value
#' @param mc.cores number of cores
#' @param mode whether to use network smoothing index (S) or network smoothing values (Xs)
#' @param ... additional parameteres of ND
#' @return \code{data.frame} with X0, Xs, S, p and Sp
#' @import parallel
#' @export
#'
calc_adjND <- function(X0, W, eps=rep(1, ncol(X0)), k=99, mode=c("S", "Xs"), mc.cores=1, return.perm=FALSE, ...){

  mode <- match.arg(mode)

  all_X0 <- c(list(X0), lapply(1:k, function(x) matrix(as.numeric(X0), ncol=ncol(X0), dimnames = list(sample(rownames(X0), nrow(X0))))))
  all_X0 <- lapply(all_X0, function(x) x[match(rownames(W), rownames(x)), , drop=F ])

  cat("network propagation\n")
  if(mc.cores==1){
    all_Xs <- lapply(all_X0, function(x) ND(x, W, ...)$Xs)
  }else{
    all_Xs <- parallel::mclapply(all_X0, function(x) ND(x, W, ...)$Xs, mc.cores=mc.cores)
  }


  if(mode=="S"){
    cat("calculation of Sp\n")

    #S
    all_S <- lapply(1:(k+1), function(i) nsi(all_X0[[i]], all_Xs[[i]], eps=eps))
    Xs <- all_Xs[[1]]

  }else{
    all_S <- all_Xs
  }

  ## trash all permutations
  if(!return.perm){
    rm(all_X0, all_Xs)
  }

  #p
  est_p <- calc_p(all_S)


  S <- all_S[[1]]

  if(mode=="S"){

    if(return.perm){

      out <- list(Xs=Xs, eps=eps, S=S, p=est_p, Sp= S * -log10(est_p), all_X0=all_X0, all_Xs=all_Xs)

    }else{

      out <- list(Xs=Xs, eps=eps, S=S, p=est_p, Sp= S * -log10(est_p))

    }

  }else{

    if(return.perm){

      out <- list(Xs=S, eps=eps, p=est_p, Xsp= S * -log10(est_p), all_X0=all_X0, all_Xs=all_Xs)

    }else{
      out <- list(Xs=S, eps=eps, p=est_p, Xsp= S * -log10(est_p))
    }
  }

  rm(all_S) ## trash all permutations

  return(out)


}
