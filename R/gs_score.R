#' Calculation of gene set scores
#' 
#' @description Calculation of gene set scores
#' @param x input vector
#' @param gsList gene set listy gene set list
#' @return numeric vector of gene set scores
#' @export
#'

gs_score <- function(x=NULL, gsList=NULL){
    return(unlist(lapply(gsList, function(xx) sum(x[xx]))))
}

