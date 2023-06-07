#' Network Smoothing Index
#' 
#' @description Network Smoothing Index
#' @param X0 input matrix/vector
#' @param Xs smoothed matrix/vector
#' @param eps numeric value
#' @return network smoothing index S
#' @usage nsi(X0, Xs, eps=rep(1, ncol(X0)))
#' @examples 
#' \dontrun{nsi(X0, Xs, eps=rep(1, ncol(X0)))}
#' @export
nsi <- function(X0=NULL, Xs=NULL, eps=NULL){

    if(!(identical(rownames(X0), rownames(Xs)))) {
        stop('rownames(X0) and rownames(Xs) are not identicail\n')     
    }
    
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

    eps <- matrix(eps, nrow=nrow(X0), ncol=ncol(X0), byrow = T)
    S <-  Xs / (X0 + eps)
    return(S)
}
