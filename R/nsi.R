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
nsi <- function(X0, Xs, eps=rep(1, ncol(X0))){

    if(!(identical(rownames(X0), rownames(Xs)))) {
        stop('rownames(X0) and rownames(Xs) are not identicail\n')     
    }

    eps <- matrix(eps, nrow=nrow(X0), ncol=ncol(X0), byrow=TRUE)
    S <-  Xs / (X0 + eps)
    return(S)
}
