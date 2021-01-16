#' assess_centrality
#'
#' @export
#' @import igraph

assess_centrality <- function(subnetwork, graph){

  ans <- data.frame(
    id=V(subnetwork)$name,
    degree=igraph::degree(subnetwork),
    betwenness=igraph::betweenness(subnetwork),
    closeness=igraph::closeness(subnetwork),
    degree_b=igraph::degree(graph, v = V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
    betwenness_b=igraph::betweenness(graph, v = V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
    closeness_b=igraph::closeness(graph, v = V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
    stringsAsFactors=F
  )


  return(ans)
}
