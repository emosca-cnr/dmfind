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

    X0_1 <- X0[, classes==1, drop=F]
    X0_2 <- X0[, classes==2, drop=F]
    Xs_1 <- Xs[, classes==1, drop=F]
    Xs_2 <- Xs[, classes==2, drop=F]

    nsi_2 <- nsi(X0_2, Xs_2, eps=eps)
    nsi_1 <- nsi(X0_1, Xs_1, eps=eps)

    delta_S <- nsi_2 - nsi_1

    dS_df <- data.frame(X0_1=as.numeric(X0_1), X0_2=as.numeric(X0_2),  
                        dX0=as.numeric(X0_2 - X0_1), Xs_1=as.numeric(Xs_1), 
                        Xs_2=as.numeric(Xs_2), dXs=as.numeric(Xs_2 - Xs_1), 
                        S_1=as.numeric(nsi_1), S_2=as.numeric(nsi_1), 
                        dS=as.numeric(delta_S), row.names = rownames(X0), 
                        stringsAsFactors = FALSE)

    return(dS_df)

}
