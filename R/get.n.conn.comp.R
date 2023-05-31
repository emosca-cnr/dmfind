#' Get the n-sized connected components
#' @param x igraph object
#' @param n network component minimum size
#' @export
#' @import igraph

get.n.conn.comp <- function(x=NULL, n=NULL){
  
  conn <- clusters(x)
  clstr.ok <- which(conn$csize >= n)
  nodes.ok <- which(conn$membership %in% clstr.ok)
  
  return(induced.subgraph(x, nodes.ok))
  
}
