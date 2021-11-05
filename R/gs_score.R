#' Calculation of gene set scores
#' 
#' @description Calculation of gene set scores
#' @param x input vector
#' @param gs_list gene set listy gene set list
#' @return numeric vector of gene set scores
#' @export
#'

gs_score <- function(x=NULL, gs_list=NULL){
    return(unlist(lapply(gs_list, function(xx) sum(x[xx]))))
}

