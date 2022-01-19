#' Extraction and plotting of connected subnetworks from a graph 
#' and a set of selected vertices
#' 
#' @description Extraction and plotting of connected subnetworks 
#' from a graph and a set of selected vertices
#' @param graph igraph object
#' @param colorBy object to be colored by an attribute
#' @param colorQuant TRUE/FALSE
#' @param labelBy select a feature to label the graph
#' @param pal palette to color the graph
#' @param community communities from network analysis 
#' @param commWin distance of the vertices inside a community
#' @param commWb distance between communities
#' @param lo layout 
#' @param legendOff legend
#' @param plotOutfile output file
#' @param ... additional arguments
#' @return plot network
#' @usage plot_network(graph, colorBy=NULL, colorQuant=TRUE, labelBy="name", 
#' pal=NULL, plotOutfile="graph.jpg", community=NULL, commWin=2, commWb=1, 
#' lo=NULL, legendOff=FALSE, ...)
#' @examples 
#' \dontrun{if (require("igraph")) plot_network(graph, colorBy=NULL, 
#' colorQuant=TRUE, labelBy="name", pal=NULL, plotOutfile="graph.jpg", 
#' community=NULL, commWin=2, commWb=1, lo=NULL, legendOff=FALSE)}
#' @export
#' @import igraph
#' @import RColorBrewer
#' @import grDevices
#' @import graphics
#' @import utils

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
        ew <- edge_weights(community, graph, weightWithin=commWin,
                           weightBetween=commWb)
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
