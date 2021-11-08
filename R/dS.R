#' Delta Network Smoothing Index
#' 
#' @description Delta Network Smoothing Index
#' @param X0 input matrix
#' @param Xs smoothed matrix
#' @param classes vector of classes, values must be equal to 1 or 2
#' @param eps numeric value
#' @return named vector of delta smoothing indexes
#'
dS <- function(X0, Xs, classes, eps=1){
    
    if(!(identical(rownames(X0), rownames(Xs))))
        stop('rownames(X0) and rownames(Xs) are not identicail\n')

    X01 <- X0[, classes==1, drop=FALSE]
    X02 <- X0[, classes==2, drop=FALSE]
    Xs1 <- Xs[, classes==1, drop=FALSE]
    Xs2 <- Xs[, classes==2, drop=FALSE]

    nsi2 <- nsi(X02, Xs2, eps=eps)
    nsi1 <- nsi(X01, Xs1, eps=eps)

    delta_S <- nsi2 - nsi1

    dSdf <- data.frame(X01=as.numeric(X01), X02=as.numeric(X02),  
                        dX0=as.numeric(X02 - X01), Xs1=as.numeric(Xs1), 
                        Xs2=as.numeric(Xs2), dXs=as.numeric(Xs2 - Xs1), 
                        S1=as.numeric(nsi1), S2=as.numeric(nsi1), 
                        dS=as.numeric(delta_S), row.names = rownames(X0), 
                        stringsAsFactors = FALSE)

    return(dSdf)

}
