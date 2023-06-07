#' Delta Network Smoothing Index
#'
#' @description Delta Network Smoothing Index
#' @param X0 input matrix
#' @param Xs smoothed matrix
#' @param classes vector of classes, values must be equal to 1 or 2
#' @param eps numeric value
#' @return named vector of delta smoothing indexes
#' @usage dS(X0, Xs, classes, eps=1)
#' @examples
#' \dontrun{dS(X0, Xs, classes, eps=1)}
#'
dS <- function(X0 = NULL,
               Xs = NULL,
               classes = NULL,
               eps = NULL) {
    
    if (!(identical(rownames(X0), rownames(Xs))))
        stop('rownames(X0) and rownames(Xs) are not identicail\n')
    
    if (is.null(eps)) {
        cat("eps value(s) not provided. It is important to tune it! Using 1...\n")
        eps <- matrix(1, nrow = 1, ncol = 2)
    }
    
    if (!is.matrix(eps)) {
        stop("eps must be a 1-row matrix with 2 columns\n")
    }
    
    if (ncol(eps) != 2) {
        stop("eps must have 2 columns\n")
    }
    
    X01 <- cbind(X01=Matrix::rowMeans(X0[, classes == 1, drop = FALSE]))
    X02 <- cbind(X02=Matrix::rowMeans(X0[, classes == 2, drop = FALSE]))
    
    Xs1 <- cbind(Xs1=Matrix::rowMeans(Xs[, classes == 1, drop = FALSE]))
    Xs2 <- cbind(Xs2=Matrix::rowMeans(Xs[, classes == 2, drop = FALSE]))
    
    nsi1 <- nsi(X01, Xs1, eps = eps[, 1, drop=F])
    nsi2 <- nsi(X02, Xs2, eps = eps[, 2, drop=F])
    
    delta_S <- nsi2 - nsi1
    
    dSdf <-
        data.frame(
            S1 = as.numeric(nsi1),
            S2 = as.numeric(nsi2),
            dS_2vs1 = as.numeric(delta_S),
            row.names = rownames(X0),
            stringsAsFactors = FALSE
        )
    
    return(dSdf)
    
}
