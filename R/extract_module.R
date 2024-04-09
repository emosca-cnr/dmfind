#' Extraction and plotting of connected subnetworks from a graph and a set of selected vertices
#' 
#' @description Extraction and plotting of connected subnetworks from a 
#' graph and a set of selected vertices
#' @param graph igraph object
#' @param NSIRes result of calc_adjND() or calc_dS()
#' @param column which column to consider among those available in the elements of calc_adjND()
#' @param selectedVertices subset of vertices
#' @param minSubnetSize minimum size of subnetwork
#' @return subnetwork 
#' @export
#' @importFrom igraph induced_subgraph V V<- clusters
#' @importFrom NPATools get_nconn_comp

extract_module <- function(graph=NULL, selectedVertices=NULL, NSIRes=NULL, column=1, minSubnetSize=2){

    #network between the selected gene
    dm <- induced_subgraph(graph, 
          V(graph)$name[match(selectedVertices, V(graph)$name)])

  
    if(any(names(NSIRes) == "Sp")){
        V(dm)$Sp <- NSIRes$Sp[match(V(dm)$name, rownames(NSIRes$Sp)), column]
        V(dm)$S <- NSIRes$S[match(V(dm)$name, rownames(NSIRes$S)), column]
        V(dm)$p <- NSIRes$p[match(V(dm)$name, rownames(NSIRes$p)), column]
        V(dm)$X0 <- NSIRes$X0[match(V(dm)$name, rownames(NSIRes$X0)), column]
    }
    
    if(any(names(NSIRes) == "dS")){
        V(dm)$dS <- NSIRes$dS[match(V(dm)$name, rownames(NSIRes$dS)), column]
        V(dm)$S1 <- NSIRes$S[match(V(dm)$name, rownames(NSIRes$S)), 1]
        V(dm)$S2 <- NSIRes$S[match(V(dm)$name, rownames(NSIRes$S)), 2]
        V(dm)$X01 <- NSIRes$X0[match(V(dm)$name, rownames(NSIRes$X0)), 1]
        V(dm)$X02 <- NSIRes$X0[match(V(dm)$name, rownames(NSIRes$X0)), 2]
    }
    

    #get components and plot the graphs
    dmClusters <- clusters(dm)
    V(dm)$subnetId <- dmClusters$membership
    V(dm)$subnetSize <- dmClusters$csize[dmClusters$membership]

    #plot
    dmComponents <- get_nconn_comp(x = dm, n = minSubnetSize)

    return(dmComponents)

}
