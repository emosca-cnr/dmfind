#' Extract the maximum connected components
#' 
#' @description Extract the maximum connected components
#' @param x network
#' @param verbose TRUE/FALSE
#' @export

get.max.conn.comp <- function(x, verbose=TRUE){
    conn <- clusters(x)
    clstr.ok <- which.max(conn$csize)
    if(verbose){
        print(table(conn$csize))
        print(head(conn$csize))
    }
    out <- induced.subgraph(x, V(x)$name[conn$membership %in% clstr.ok])
    return(out)
}
