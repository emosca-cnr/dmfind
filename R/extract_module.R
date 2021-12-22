#' Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' 
#' @description Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' @param graph igraph object
#' @param selectedVertices subset of vertices
#' @param minSubnetSize minimum size of subnetwork
#' @return subnetwork 
#' @usage extract_module(graph, selectedVertices, X0=NULL, minSubnetSize=2)
#' @examples extract_module(graph, selectedVertices, X0=NULL, minSubnetSize=2)
#' @export
#' @import igraph

extract_module <- function(graph, selectedVertices, 
                           X0=NULL, minSubnetSize=2){

    #network between the selected gene
    dm <- igraph::induced_subgraph(graph, 
          V(graph)$name[match(selectedVertices, V(graph)$name)])

    if(!is.null(X0)){
        V(dm)$X0 <- X0[match(V(dm)$name, rownames(X0)), 1]
    }

    #get components and plot the graphs
    dmClusters <- igraph::clusters(dm)
    V(dm)$subnetId <- dmClusters$membership
    V(dm)$subnetSize <- dmClusters$csize[dmClusters$membership]

    #plot
    dmComponents <- igraph::induced_subgraph(dm, 
                     V(dm)$name[V(dm)$subnetSize >= minSubnetSize])

    return(dmComponents)

}
