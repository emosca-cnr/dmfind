#' Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' 
#' @description Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' @param graph igraph object
#' @param selected_vertices subset of vertices
#' @param min_subnet_size minimum size of subnetwork
#' @export
#' @import igraph

extract_module <- function(graph, selected_vertices, 
                           X0=NULL, min_subnet_size=2){

    #network between the selected gene
    dm <- igraph::induced_subgraph(graph, 
          V(graph)$name[match(selected_vertices, V(graph)$name)])

    if(!is.null(X0)){
        V(dm)$X0 <- X0[match(V(dm)$name, rownames(X0)), 1]
    }

    #get components and plot the graphs
    dm_clusters <- igraph::clusters(dm)
    V(dm)$subnet_id <- dm_clusters$membership
    V(dm)$subnet_size <- dm_clusters$csize[dm_clusters$membership]

    #plot
    dm_components <- igraph::induced_subgraph(dm, 
                     V(dm)$name[V(dm)$subnet_size >= min_subnet_size])

    return(dm_components)

}
