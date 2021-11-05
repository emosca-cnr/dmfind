#' Calculation of delta network smoothing indexes
#' 
#' @description Calculation of delta network smoothing indexes
#' @param X0 input matrix
#' @param A adjancency matrix
#' @param classes numeric vector of \{1,2\}
#' @param eps numeric value
#' @param ... additional parameteres of ND
#' @return \code{data.frame} with delta gene network smoothing indexes
#' @export
#'
calc_dS <- function(X0, A, classes, eps=1, ...){

  cat("network propagation\n")
  Xs <- ND(X0, A, ...)$Xs

  cat("calculation of dS\n")
  dS <- vector('list', 2)

  dS_df <- dS(X0, Xs, classes=classes, eps=eps)

  return(dS_df)

}
