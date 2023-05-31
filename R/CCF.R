CCF <- function(x){
	
	ans <- igraph::clusters(x)
	ans <- sum(ans$csize[ans$csize>1]) / length(ans$membership)
	
	return(ans)
	
}