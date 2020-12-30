#' Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' @param graph igraph object
#' @param selected_vertices subset of vertices
#' @param vertices_subset a subset of selected_vertices (optional)
#' @param vertices_weight weights (optional)
#' @param vertices_label labels (optional)
#' @param min_subnet_size minimum size of subnetwork
#' @param plot_outfile output file
#' @param plot_flag whether to plot or not the network
#' @param ... additional arguments
#' @export
#' @import igraph
#' @import RColorBrewer

extract_module <- function(graph, selected_vertices, vertices_subset=NULL, vertices_weight=NULL, vertices_label=NULL, min_subnet_size = 2, plot_flag=TRUE, plot_outfile='graph.jpg', ...){

  #network between the selected gene
  dm <- induced_subgraph(graph, V(graph)$name[match(selected_vertices, V(graph)$name)])

  #get components and plot the graphs
  V(dm)$subnet_id <- clusters(dm)$membership
  V(dm)$subnet_size <- clusters(dm)$csize[clusters(dm)$membership]

  #plot
  dm_components <- induced_subgraph(dm, V(dm)$name[V(dm)$subnet_size >= min_subnet_size])

  #shape
  V(dm_components)$shape <- rep('circle', length(V(dm_components)))
  if(length(vertices_subset) > 0){
    V(dm_components)$shape[V(dm_components)$name %in% vertices_subset] <- 'square'
  }
  #color
  if(length(vertices_weight) >= length(V(dm_components))){
    idx <- match(V(dm_components)$name, names(vertices_weight))
    if(any(is.na(idx)))
      stop('vertices weights do not exist for every vertex')
    color_values <- cut(vertices_weight[match(V(dm_components)$name, names(vertices_weight))], 9)
    V(dm_components)$color <- RColorBrewer::brewer.pal(9, "Purples")[as.numeric(color_values)]
  }
  #label
  if(length(vertices_label) >0){
    V(dm_components)$label <- vertices_label[match(V(dm_components)$name, names(vertices_label))]
  }else{
    V(dm_components)$label <- V(dm)$subnet_id[match(V(dm_components)$name, V(dm)$name)]
  }

  #same layout for both figures
  if(plot_flag==TRUE){
    lo <- layout.fruchterman.reingold(dm_components)

    jpeg(plot_outfile, width = 200, height = 200, units='mm', res=300)
    par(mar=c(1, 1, 1, 1))
    plot.igraph(dm_components, layout=lo, ...)
    if(length(vertices_weight) >= length(V(dm_components))){
      legend("bottomright", levels(factor(color_values, levels=sort(unique(color_values)))), col=RColorBrewer::brewer.pal(9, "Purples"), pch=16, bty="n", cex=0.6)
    }
    dev.off()
  }

  return(list(full=dm, conn_subnets=dm_components))

}
