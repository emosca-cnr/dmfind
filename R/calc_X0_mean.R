#' Calculate X0 mean
#' @param X0 vector or matrix composed of column vectors with initial 
#' distribution of information
#' @param classes vector of class labels for the columns of X0. Values must be
#' equal to 1 or 2
#' 
#' @export
#' 
calc_X0_mean <- function(X0=NULL, classes=NULL){
  
  X0_means <- t(apply(X0, 1, function(i_row) tapply(i_row, INDEX = classes, FUN = mean)))
  
  return(X0_means)
  
}