#'
#'
#'
#'
#'

assess_centrality <- function(subnetwork, graph){


  ans <- data.frame(
    id=V(subnetwork)$name,
    degree=degree(subnetwork),
    betwenness=betweenness(subnetwork),
    closeness=closeness(subnetwork),
    degree_b=degree(graph, v = V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
    betwenness_b=betweenness(graph, v = V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
    closeness_b=closeness(graph, v = V(graph)$name[match(V(subnetwork)$name, V(graph)$name)]),
    stringsAsFactors=F
  )


  return(ans)
}
