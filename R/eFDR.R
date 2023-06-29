#' empirical False Discovery Rate
#' @param realValues vactor of real (observed) values
#' @param realValues vactor of all values (real + permuted) values
#' @param mc.cores number of cores
#' @import parallel
eFDR <- function(realValues, allValues, BPPARAM = NULL){
  
  #FDR = (# x > permuted values / # x > real values) * (#real values / #permuted values)
  #approch used in GSEA
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  fdr_values <- length(realValues) / length(allValues)
  fdr_values <- fdr_values * unlist(bplapply(realValues, function(x) sum(allValues >= x) / sum(realValues >= x), BPPARAM = BPPARAM))
  
  #bound to 1
  fdr_values[fdr_values>1] <- 1 

  return(fdr_values)
  
}