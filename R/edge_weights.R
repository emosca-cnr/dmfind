#' Function for the assignment of edges weights in order to better visualize the communities
#' 
#' @description Function for the assignment of edges weights in order to better visualize the communities
#' @param community list with the mandatory field "membership"
#' @param network pathway network
#' @param weightWithin value to weight the attraction between two vertices of the same community
#' @param weightBetween value to weight the attraction between two vertices of two disntict communities
#' @return vector of edge weights
#' @import igraph
#' @export
#'
edge_weights <- function(community, network, weightWithin=100, 
                         weightBetween=1) {
    bridges <- igraph::crossing(communities=community, graph=network)
    weights <- ifelse(test=bridges, yes=weightBetween, no=weightWithin)
    return(weights)
}
