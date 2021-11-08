#' assess centrality measure to the top network
#' 
#' @param subnetwork portion of graph, top network
#' @param graph igraph object, total graph
#' @description Function for assess centrality measures to the top network.
#' @usage assess_centrality(subnetwork, graph)
#' @example assess_centrality(subnetwork, graph)
#' @export
#' @import igraph

assess_centrality <- function(subnetwork, graph){

    ans <- data.frame(id=V(subnetwork)$name,
            degree=igraph::degree(subnetwork),
            betwenness=igraph::betweenness(subnetwork),
            closeness=igraph::closeness(subnetwork),
            degree_b=igraph::degree(graph, 
            v=V(graph)$name[match(V(subnetwork)$name,V(graph)$name)]),
            betwenness_b=igraph::betweenness(graph, 
            v=V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
            closeness_b=igraph::closeness(graph,
            v=V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
            stringsAsFactors=F)

    return(ans)
}
