#' Calculation of pathway score matrix
#' @param x input matrix
#' @param gs_list gene set list
#' @param use_median whether to use median (TRUE) or mean (FALSE)
#' @param scaling_eps numeric value
#' @param mc.cores number of cores
#' @return matrix of scaled pathway scores
#' @export
#'
std_gs_score <- function(x, gs_list, use_median=FALSE, scaling_eps=1, mc.cores=1){

  #calculation of gs_score and scaled gs_score
  out <- split(x, rep(1:ncol(x), each = nrow(x)))

  if(mc.cores > 1){
    out <- as.data.frame(parallel::mclapply(out, gs_score, gs_list=gs_list, mc.cores=mc.cores))
  }else{
    out <- as.data.frame(lapply(out, gs_score, gs_list=gs_list))
  }

  colnames(out) <- colnames(x)

  if(use_median){
    out <- t(scale(t(out), center=apply(out, 1, median), scale=apply(out, 1, sd) + scaling_eps))
  }else{
    out <- t(scale(t(out), center=apply(out, 1, mean), scale=apply(out, 1, sd) + scaling_eps))
  }

  return(out)

}
