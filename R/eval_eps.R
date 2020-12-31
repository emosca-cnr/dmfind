#' Evaluate the effect of epsilon on NSI
#' @param X0 initial condition matrix
#' @param Xs network smoothing result
#' @param eps matrix with epsilon values (rows) for each column of X0, Xs
#' @param G matrix with epsilon values (rows) for each column of X0, Xs
#' @param top matrix with epsilon values (rows) for each column of X0, Xs
#' @return list with sn_g and sn_summary elements
#' @import Bolstad2
#' @export
#'

eval_eps <- function(X0, Xs, eps=NULL, G, top=min(300, nrow(X0))){

  if(is.null(eps)){
    stop("eps is mandatory")
  }


  nsi_eps <- vector("list", nrow(eps))
  csX0 <- nsi_eps
  csXs <- nsi_eps
  for(i in 1:nrow(eps)){

    temp_eps <- matrix(eps[i, ], ncol=ncol(X0), nrow = nrow(X0), byrow = T)

    #NSI
    nsi_eps[[i]] <- Xs / (X0 + temp_eps)

    #order by current NSI
    nsi_order <- apply(nsi_eps[[i]], 2, function(x) order(-x)[1:top])

    #omega on X0
    csX0[[i]] <- sapply(1:ncol(X0), function(j) omega(G, array(X0[, j], dimnames = list(rownames(X0)))[nsi_order[, j]]))

    #omega on Xs
    csXs[[i]] <- sapply(1:ncol(Xs), function(j) omega(G, array(Xs[, j], dimnames = list(rownames(Xs)))[nsi_order[, j]]))

  }

  aucX0 <- unlist(lapply(out$csX0, function(x) Bolstad2::sintegral(x = 1:top, fx = x)$int))
  aucXs <- unlist(lapply(out$csXs, function(x) Bolstad2::sintegral(x = 1:top, fx = x)$int))

  out <- list(S=nsi_eps, csX0=csX0, csXs=csXs, aucX0=aucX0, aucXs=aucXs)

  return(out)

}
