#' Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' 
#' @description Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' @param graph igraph object
#' @param selectedVertices subset of vertices
#' @param verticesSubset a subset of selectedVertices (optional)
#' @param verticesWeight weights (optional)
#' @param verticesLabel labels (optional)
#' @param minSubnet_Size minimum size of subnetwork
#' @param plotOutfile output file
#' @param plotFlag whether to plot or not the network
#' @param ... additional arguments
#' @export
#' @import igraph
#' @import RColorBrewer

plot_network <- function(graph, colorBy=NULL, colorQuant=TRUE, 
                         labelBy="name", pal=NULL, plotOutfile="graph.jpg", 
                         community=NULL, commWin=2, commWb=1, lo=NULL, 
                         legendOff=FALSE, ...){

    if(!is.null(colorBy)){
        colorValues <- get.vertex.attribute(graph, colorBy)
        if(colorQuant){
            colorValues <- cut(colorValues, length(pal), dig.lab = 1)
        }
        V(graph)$color <- pal[as.numeric(as.factor(colorValues))]
    }
    if(labelBy == ""){
        vertexLabels <- ""
    }else{
        vertexLabels <- get.vertex.attribute(graph, labelBy)
     }
    V(graph)$label <- vertexLabels

    ew <- NULL
    if(!is.null(community)){
        ew <- edge_weights(community, graph, weight.within=commWin,
                           weight.between=commWb)
    }

    if(is.null(lo)){
        lo <- layout_with_fr(graph, weights=ew)
    }
  
    if(!is.null(plotOutfile)){
        jpeg(plotOutfile, width=210, height = 200, units='mm', res=300)
    }
  
    if(!is.null(colorBy) & !legendOff){
        layout(matrix(c(1, 2), ncol=2), widths = c(0.9, 0.2))
    }
    par(mar=c(1, 1, 1, 1))

    plot.igraph(graph, vertex.labels=vertexLabels, layout=lo, ...)

    if(!is.null(colorBy) & !legendOff){
        plot.new()
        legend("bottomright", levels(factor(colorValues, 
              levels=sort(unique(colorValues)))), 
              col=pal, pch=16, bty="n", cex=0.6)
    }

    if(!is.null(plotOutfile)){
        dev.off()
    }
  
    return(lo)

}
