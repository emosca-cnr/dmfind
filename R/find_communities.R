#' Find topological communities
#' 
#' @description Find topological communities
#' @param g graph
#' @param eWeights edge weights
#' @param vWeights vertex weights
#' @param verbose TRUE/FALSE
#' @param methods one or more of "fastgreedy", "labprop", "walktrap", 
#' "eigen", "multilev", "infomap"
#' @return list od community objects
#' @usage find_communities(g, eWeights=NULL, vWeights=NULL, verbose=TRUE, 
#' methods=c("fastgreedy", "labprop", "walktrap", "eigen", 
#' "multilev", "infomap"))
#' @examples  
#' \dontrun{find_communities(g, eWeights=NULL, vWeights=NULL, verbose=TRUE,
#'  methods=c("fastgreedy", "labprop", "walktrap", "eigen",
#'   "multilev", "infomap"))}
#' @import igraph
#' @export
find_communities <- function(g=NULL, eWeights=NULL, vWeights=NULL, 
                             verbose=TRUE, methods=c("fastgreedy", "multilev")){

    commList <- setNames(vector('list', length = length(methods)), methods)
    
    commInfo <- data.frame(algorithm=methods,
                            modularity=NA,
                            n=NA,
                            stringsAsFactors=FALSE)

    if("fastgreedy" %in% methods){
        if(verbose)
            print("fastgreedy")
        idx <- which(grepl("fastgreedy", methods))
        commList[[idx]] <- igraph::fastgreedy.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            igraph::modularity(g, commList[[idx]]$membership,
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("labprop" %in% methods){
        if(verbose)
            print("labprop")
        idx <- which(grepl("labprop", methods))
        commList[[idx]] <- 
            igraph::label.propagation.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            igraph::modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("walktrap" %in% methods){
        if(verbose)
            print("walktrap")
        idx <- which(grepl("walktrap", methods))
        commList[[idx]] <- igraph::walktrap.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            igraph::modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("eigen" %in% methods){
        if(verbose)
            print("eigen")
        idx <- which(grepl("eigen", methods))
        commList[[idx]] <- 
            igraph::leading.eigenvector.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            igraph::modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
    if("multilev" %in% methods){
        if(verbose)
            print('multilev')
        idx <- which(grepl("multilev", methods))
        commList[[idx]] <-
            igraph::multilevel.community(g, weights=eWeights)
        commInfo$modularity[idx] <- 
            igraph::modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
    }
  if("infomap" %in% methods){
        if(verbose)
            print('infomap')
        idx <- which(grepl("infomap", methods))
        commList[[idx]] <- 
            igraph::infomap.community(g, eWeights=eWeights, 
                                      vWeights=vWeights)
        commInfo$modularity[idx] <- 
            igraph::modularity(g, commList[[idx]]$membership, 
                               weights=eWeights)
        commInfo$n[idx] <- max(commList[[idx]]$membership)
  }

  return(list(comm=commList, info=commInfo))

}
