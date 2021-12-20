#' Extract a connected components of size n
#' 
#' @description Extract the maximum connected components
#' @param x network
#' @param n size of connected component 
#' @param verbose TRUE/FALSE
#' @return subgraph
#' @usage get_n_conn_comp(x,n)
#' @examples get_n_conn_comp(x,n)
#' @export

get_n_conn_comp <- function(x, n){
    conn <- clusters(x)
    clstrOk <- which(conn$csize >= n)
    nodesOk <- which(conn$membership %in% clstrOk)
    return(induced.subgraph(x, nodesOk))
}
