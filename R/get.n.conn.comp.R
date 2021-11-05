#' Extract a connected components of size n
#' 
#' @description Extract the maximum connected components
#' @param x network
#' @param n size of connected component 
#' @param verbose TRUE/FALSE
#' @export

get.n.conn.comp <- function(x, n){
    conn <- clusters(x)
    clstr.ok <- which(conn$csize >= n)
    nodes.ok <- which(conn$membership %in% clstr.ok)
    return(induced.subgraph(x, nodes.ok))
}
