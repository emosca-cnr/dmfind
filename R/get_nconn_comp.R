#' Extract a connected components of size n
#' 
#' @description Extract the maximum connected components
#' @param x network
#' @param n size of connected component 
#' @return subgraph
#' @usage get_nconn_comp(x,n)
#' @examples 
#' \dontrun{get_nconn_comp(x,n)}
#' @import igraph
#' @export

get_nconn_comp <- function(x, n){
    conn <- clusters(x)
    clstrOk <- which(conn$csize >= n)
    nodesOk <- which(conn$membership %in% clstrOk)
    return(induced.subgraph(x, nodesOk))
}
