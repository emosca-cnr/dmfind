#' Connected component fraction
#' @export
#' @param x igraph object
#' @param min_size min size of a connected component
#' @importFrom igraph clusters
CCF <- function(x=NULL, min_size=2){
	
	ans <- clusters(x)
	ans <- sum(ans$csize[ans$csize>1]) / length(ans$membership)
	
	return(ans)
	
}