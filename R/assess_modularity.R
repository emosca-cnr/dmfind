#' Assess modularity
#' @export
#' @import BiocParallel
#' 
assess_modularity <- function(G = NULL, vertices = NULL, ranks = NULL, commMethods=c("fastgreedy", "multilev"),  BPPARAM = NULL) {
  
  if (is.null(ranks)) {
    ranks <- seq(10, length(vertices), by = 10)
  }
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  
  top_nets <- lapply(ranks, function(k) induced_subgraph(graph = G, vids = V(G)$name[V(G)$name %in% vertices[1:k]]))
  names(top_nets) <- ranks
  
  top_comm <- BiocParallel::bplapply(top_nets, function(g_i) {
    #print(x)
    comm <- find_communities(g_i, verbose = F, methods = commMethods)
    best <- comm$info
    best <- best[order(-best$modularity, best$n),]
    tmp <- array(c(vcount(g_i), ecount(g_i), best[1, ]))
    return(tmp)
  }, BPPARAM = BPPARAM)
  
  top_comm <- do.call(rbind, top_comm)
  df.out <- data.frame(rank = names(top_nets),
                       n_vertex = as.numeric(top_comm[, 1]),
                       n_edge = as.numeric(top_comm[, 2]),
                       algorithm = unlist(top_comm[, 3]),
                       modularity = as.numeric(top_comm[, 4]),
                       n_community = as.numeric(top_comm[, 5]), stringsAsFactors = F)
  
  
  return(df.out)
}