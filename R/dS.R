#' Delta Network Smoothing Index
#'
#' @description Delta Network Smoothing Index
#' @param X0 input matrix
#' @param Xs smoothed matrix
#' @param classes vector of classes, values must be equal to 1 or 2
#' @param eps numeric value
#' @return named vector of delta smoothing indexes
#' @importFrom Matrix Matrix rowMeans
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
    
    classes <- as.numeric(as.factor(classes))
    
    X01 <- Matrix(rowMeans(X0[, classes == 1, drop = FALSE]), sparse = T, dimnames = list(rownames(X0), "X01"))
    X02 <- Matrix(rowMeans(X0[, classes == 2, drop = FALSE]), sparse = T, dimnames = list(rownames(X0), "X02"))
    
    Xs1 <- Matrix(rowMeans(Xs[, classes == 1, drop = FALSE]), sparse = T, dimnames = list(rownames(Xs), "Xs1"))
    Xs2 <- Matrix(rowMeans(Xs[, classes == 2, drop = FALSE]), sparse = T, dimnames = list(rownames(X0), "Xs2"))
    
    nsi1 <- nsi(X01, Xs1, eps = eps[, 1, drop = F])
    nsi2 <- nsi(X02, Xs2, eps = eps[, 2, drop = F])
    
    delta_S <- nsi2 - nsi1
    
    ans <- list(
        X0 = cbind(X01, X02),
        Xs = cbind(Xs1, Xs2),
        eps = eps,
        S = data.frame(S1=as.numeric(nsi1), S2=as.numeric(nsi1), row.names = rownames(nsi1)),
        dS = data.frame(dS=as.numeric(delta_S), row.names = rownames(delta_S))
    )
    
    return(ans)
    
}
