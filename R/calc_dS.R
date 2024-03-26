#' Calculation of delta network smoothing indexes
#' 
#' @description Calculation of delta network smoothing indexes
#' @param X0 input matrix
#' @param A adjancency matrix
#' @param classes numeric vector of \{1,2\}
#' @param eps 1-row matrix with two columns for the eps values of class 1 and 2
#' respectively
#' @param ... additional parameteres of ND
#' @usage calc_dS(X0, A, classes, eps=1, ...)
#' @examples 
#' \dontrun{calc_dS(X0, A, classes, eps=1)}
#' @return \code{data.frame} with delta gene network smoothing indexes
#' @export
#' @import NPATools
calc_dS <- function(X0=NULL, W=NULL, classes=NULL, eps=NULL, ...){

    if(is.null(classes)){
        stop("missing mandatory input 'classes'\n")
    }
    if (is.null(eps)) {
        cat("eps value(s) not provided. It is important to tune it! Using 1...\n")
        eps <- matrix(1, nrow = 1, ncol = 2)
    }
    
    if (!is.matrix(eps)) {
        stop("eps must be a 1-row matrix with two columns\n")
    }
    
    if (ncol(eps) != 2) {
        stop("eps must have 2 columns\n")
    }
    
    cat("network propagation\n")
    Xs <- ND(X0 = X0, W = W, ...)

    cat("calculation of dS\n")

    dS_df <- dS(X0, Xs, classes=classes, eps=eps)

    return(dS_df)

}
