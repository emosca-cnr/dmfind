#' Network comparison
#'
#' @param NetList interactomes, list of igraph objects
#' @param union Perform the union of the input set of networks
#' @param intersection Perform the insersection of the input set of networks
#' @param centrMeas Calculate the centrality measures for every input network
#' @param ji Calculation of the Jaccard index
#' @param oc Calculation of the overlap Coefficient
#' @return A list containing the occurrences matrix and, if selected, the union of the networks,
#' the intersection of the networks, a dataframe containing the centrality
#' measure for every network and the Jaccard Index.
#' @description Network comparison
#' @usage NetComp(NetList, union=FALSE, intersection=FALSE,
#' centrMeas=FALSE, ji=FALSE)
#' @examples
#' \dontrun{NetComp(NetList, union=FALSE, intersection=FALSE,
#' centrMeas=FALSE, ji=FALSE)}
#' @import igraph
#' @export
#'

NetComp <- function(NetList, union=FALSE, intersection=FALSE,
                    centrMeas=FALSE, ji=FALSE, oc=FALSE) {
    ### 1) create occurrences matrix (DEFAULT)
    genes <- data.frame(V(NetList[[1]])$name)
    for (i in 2:length(NetList)) {
        genes <- merge(genes, as.data.frame(V(NetList[[i]])$name), by=1, all=T)
        names(genes) <- "genes"
    }
    cat(paste0("Creating occurrences matrix for ",
               length(NetList), " interactomes...\n"))
    mOcc <- matrix(0, nrow=nrow(genes), ncol=length(NetList))
    rownames(mOcc) <- genes$genes
    mOcc <- as.data.frame(mOcc)
    for (i in 1:ncol(mOcc)) {
        mOcc[,i] <- ifelse(rownames(mOcc) %in%
                               V(NetList[[i]])$name, 1, 0)
    }
    cat("Done\n")

    ### 2) extract the union of the input interactomes
    if (union == TRUE) {
        cat(paste0("Performing union of ",
                   length(NetList), " interactomes...\n"))
        for (i in 1:length(NetList)) {
            NetList[[i]] <- NetList[[i]] %>%
                igraph::set_vertex_attr(paste0("w", i), value=1)
            NetList[[i]] <- NetList[[i]] %>%
                igraph::set_edge_attr(paste0("w", i), value=1)
        }
        if (length(NetList)>2) {
            unionNet <- igraph::union(NetList[[1]],
                                      NetList[[2]], byname=TRUE)
            for (i in 3:length(NetList)) {
                unionNet <- igraph::union(unionNet,
                                          NetList[[i]], byname=TRUE)
            }
        } else {
            unionNet <- igraph::union(NetList[[1]], NetList[[2]], byname=TRUE)
        }
        # set w as the sum of w1, ..., wn
        unionNet <- unionNet %>% igraph::set_vertex_attr("w", value=1)
        unionNet <- unionNet %>% igraph::set_edge_attr("w", value=1)
        # set missing w_i to 0
        for (i in 1:length(NetList)) {
            vertex_attr(unionNet, paste0("w",i))[is.na(igraph::vertex_attr(
                unionNet, paste0("w",i)))] <- 0
            edge_attr(unionNet, paste0("w",i))[is.na(igraph::edge_attr(
                unionNet, paste0("w",i)))] <- 0
        }
        # count the number of interactomes in which the connection exists
        if (length(NetList)>2) {
            V(unionNet)$w <- igraph::vertex_attr(unionNet, paste0("w", 1)) +
                igraph::vertex_attr(unionNet, paste0("w", 2))
            E(unionNet)$w <- igraph::edge_attr(unionNet, paste0("w", 1)) +
              igraph::edge_attr(unionNet, paste0("w", 2))
            for (i in 3:length(NetList)) {
                V(unionNet)$w <- igraph::vertex_attr(unionNet, "w") +
                    igraph::vertex_attr(unionNet, paste0("w", i))
                E(unionNet)$w <- igraph::edge_attr(unionNet, "w") +
                  igraph::edge_attr(unionNet, paste0("w", i))
            }
        } else {
            V(unionNet)$w <- igraph::vertex_attr(unionNet, paste0("w", 1)) +
                igraph::vertex_attr(unionNet, paste0("w", 2))
            E(unionNet)$w <- igraph::edge_attr(unionNet, paste0("w", 1)) +
              igraph::edge_attr(unionNet, paste0("w", 2))
        }
        table(V(unionNet)$w)
        table(E(unionNet)$w)
        cat("Done\n")
    } else {
        unionNet <- "NA"
    }

    ### 3) extract the intersection of the input interactomes
    if (intersection == TRUE) {
        cat(paste0("Extracting the intersection for ",
                   length(NetList), " interactomes...\n"))
        topEdges <-
          subgraph.edges(unionNet,
                         E(unionNet)[(E(unionNet)$w==max(E(unionNet)$w))],
                         delete.vertices = TRUE)
        topVertices <-
          induced_subgraph(unionNet,
                           V(unionNet)[(V(unionNet)$w==max(V(unionNet)$w))])
        cat("Done\n")
    } else {
        topEdges <- "NA"
    }

    ### 4) centrality measures variations
    if (centrMeas == TRUE) {
        cat(paste0("Calculating centrality measures for ",
                   length(NetList), " interactomes...\n"))
        listCentr <- list()
        for (i in 1:length(NetList)) {
            df_all <- data.frame(id=V(NetList[[i]])$name,
                      symbol=V(NetList[[i]])$label,
                      degree=igraph::degree(NetList[[i]]),
                      betwenness=igraph::betweenness(NetList[[i]]),
                      closeness=igraph::closeness(NetList[[i]]),
                      eigenCentr=igraph::eigen_centrality(NetList[[i]])$vector,
                      stringsAsFactors=FALSE)
            listCentr[[i]] <- df_all
        }
        cat("Done\n")
    } else {
        listCentr <- "NA"
    }

    ### 5) calculate Jaccard index
    if (ji == TRUE) {
        cat(paste0("Calculating Jaccard Index...\n"))
        Jindex <- length(V(topVertices))/length(V(unionNet))
        cat("Done\n")
    } else {
        Jindex <- "NA"
    }

    ### 6) calculate Overlap coefficient
    if (oc == TRUE) {
      cat(paste0("Calculating Overlap Coefficient...\n"))
      for (i in 1:length(NetList)) {
        minOv <- min(length(V(NetList[[i]])))
      }
      Overlap <- length(V(topVertices))/minOv
      cat("Done\n")
    } else {
      Overlap <- "NA"
    }

    ans <- list(mOcc, unionNet, topEdges, topVertices,
                listCentr, Jindex, Overlap)
    names(ans) <- c("Occurrence Matrix", "Network Union",
                    "Network Intersection - Edges", "Network Intersection
                    - Vertices", "Centrality Measures", "Jaccard Index",
                    "Overlap Coefficient")
    return(ans)

}

