#' reduce the network as a community network
#'
#' @param g top network, with vertex attribute "comm_id"
#' @description Function for reduce the network as a community network
#' @usage comm_net(g)
#' @examples
#' \dontrun{comm_net(g)}
#' @return dataframe with number of shared links between couple of communities
#' @export
#' @import igraph
#' 
comm_net <- function(g) {
    commNet <- igraph::contract.vertices(g, V(g)$comm_id, vertex.attr.comb = "first")
    V(commNet)$name <- V(commNet)$comm_id
    V(commNet)$label <- V(commNet)$comm_id
    
    commNet <- igraph::simplify(commNet, remove.multiple = F, remove.loops = T)
    
    E(commNet)$weight <- 1
    commNetRes <- simplify(commNet, remove.loops=FALSE)
    
    commNetDf <- data.frame(as_edgelist(commNetRes, names = TRUE),  E(commNetRes)$weight)
    colnames(commNetDf) <- c("Node 1", "Node 2", "Connections")
    l <- list(commNetRes, commNetDf)
    return(l)
    #return(commNetDf)
}
