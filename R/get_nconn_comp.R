#' Extract connected components of at least n vertices
#' 
#' @description Extract connected components of at least n vertices
#' @param x igraph obejct
#' @param n minimum size of a connected component
#' @return igraph object
#' @usage get_nconn_comp(x, n)
#' @examples 
#' \dontrun{get_nconn_comp(x,n)}
#' @import igraph
#' @export

get_nconn_comp <- function(x=NULL, n=NULL){
    conn <- clusters(x)
    clstrOk <- which(conn$csize >= n)
    nodesOk <- which(conn$membership %in% clstrOk)
    return(induced.subgraph(x, nodesOk))
}
