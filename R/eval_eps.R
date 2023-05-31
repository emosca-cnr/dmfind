#' Evaluate the effect of epsilon on NSI
#'
#' @description Evaluate the effect of epsilon on NSI
#' @param X0 initial condition matrix
#' @param Xs network smoothing result
#' @param eps matrix with epsilon values (rows) for each column of X0, Xs
#' @param G igraph object
#' @param top top genes to consider
#' @return list with sn_g and sn_summary elements
#' @usage eval_eps(X0, Xs, eps=NULL, G, top)
#' @examples
#' \dontrun{eval_eps(X0, Xs, eps=NULL, G, top)}
#' @import Bolstad2
#' @export
#'

eval_eps <- function(X0, Xs, eps=NULL, G, top=min(300, nrow(X0))){

    if(is.null(eps)){
        stop("eps is mandatory")
    }

    nsiEps <- vector("list", nrow(eps))
    csX0 <- nsiEps
    csXs <- nsiEps
    for(i in 1:nrow(eps)){
        tempEps <- matrix(eps[i, ], ncol=ncol(X0), nrow=nrow(X0), byrow=TRUE)
        #NSI
        nsiEps[[i]] <- Xs / (X0 + tempEps)
        #order by current NSI
        nsiOrder <- apply(nsiEps[[i]], 2, function(x) order(-x)[1:top])
        #omega on X0
        csX0[[i]] <- sapply(1:ncol(X0), function(j)
            omega(G, array(X0[, j],
            dimnames=list(rownames(X0)))[nsiOrder[, j]]))
        #omega on Xs
        csXs[[i]] <- sapply(1:ncol(Xs), function(j)
            omega(G, array(Xs[, j],
            dimnames=list(rownames(Xs)))[nsiOrder[, j]]))
    }

    aucX0 <- unlist(lapply(csX0, function(x)
        Bolstad2::sintegral(x=1:top, fx=x)$int))
    aucXs <- unlist(lapply(csXs, function(x)
        Bolstad2::sintegral(x=1:top, fx=x)$int))
    out <- list(S=nsiEps, csX0=csX0, csXs=csXs, aucX0=aucX0, aucXs=aucXs, eps=eps)

    return(out)

}
