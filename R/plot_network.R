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

plot_network <- function(graph, color_by=NULL, color_quant=TRUE, label_by="name", pal=NULL, plot_outfile='graph.jpg', community=NULL, comm_w_in=100, comm_w_b=1, ...){


  color_values <- get.vertex.attribute(graph, color_by)
  if(color_quant){
    color_values <- cut(color_values, length(pal))
  }
  V(graph)$color <- pal[as.numeric(color_values)]

  if(label_by == ""){
    vertex_labels <- ""
  }else{
    vertex_labels <- get.vertex.attribute(graph, label_by)
  }
  V(graph)$label <- vertex_labels

  ew <- NULL
  if(!is.null(community)){
    ew <- edge_weights(community, graph, weight.within = comm_w_in, weight.between = comm_w_b)
  }


  lo <- layout_with_fr(graph, weights = ew)

  jpeg(plot_outfile, width = 200, height = 200, units='mm', res=300)
  par(mar=c(1, 1, 1, 1))

  plot.igraph(graph, vertex.labels=vertex_labels, layout=lo, ...)

  legend("bottomright", levels(factor(color_values, levels=sort(unique(color_values)))), col=pal, pch=16, bty="n", cex=0.6)

  dev.off()

}
