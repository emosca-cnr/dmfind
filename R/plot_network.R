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

plot_network <- function(graph, color_by=NULL, color_quant=TRUE, label_by="name", pal=NULL, plot_outfile="graph.jpg", community=NULL, comm_w_in=2, comm_w_b=1, lo=NULL, legend.off=FALSE, ...){


  if(!is.null(color_by)){
    color_values <- get.vertex.attribute(graph, color_by)
    if(color_quant){
      color_values <- cut(color_values, length(pal), dig.lab = 1)
    }
    V(graph)$color <- pal[as.numeric(as.factor(color_values))]
  }
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

  if(is.null(lo)){
    lo <- layout_with_fr(graph, weights = ew)
  }
  
  if(!is.null(plot_outfile)){
    jpeg(plot_outfile, width = 210, height = 200, units='mm', res=300)
  }
  
  if(!is.null(color_by) & !legend.off){
    layout(matrix(c(1, 2), ncol = 2), widths = c(0.9, 0.2))
  }
  par(mar=c(1, 1, 1, 1))

  plot.igraph(graph, vertex.labels=vertex_labels, layout=lo, ...)

  if(!is.null(color_by) & !legend.off){
    plot.new()
    legend("bottomright", levels(factor(color_values, levels=sort(unique(color_values)))), col=pal, pch=16, bty="n", cex=0.6)
  }

  if(!is.null(plot_outfile)){
    dev.off()
  }
  

  return(lo)

}
