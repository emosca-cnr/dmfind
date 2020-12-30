#' Calculation of permutation-adjusted network smoothing index
#' @param G igraph graph object
#' @param ranked_list ranked vector; names must match V(G)$name
#' @param idx vector of ranks according to which networks will be extracted
#' @return list with sn_g and sn_summary elements
#' @import igraph
#' @export
#'

get_top_subnetworks <- function(G, ranked_list, idx=1:length(ranked_list)){


  if(is.null(idx)){
    idx_curr <- 1:length(ranked_list)
  }else{
    idx_curr <- idx
  }
  g_out <- vector('list', length = length(idx_curr))
  stat_out <- data.frame(idx=idx_curr, n_v=idx_curr, n_e=0, d=0, n_comp=0, n_v_max=0, n_v_min=0, n_v_conn_tot=0)
  j <- 1
  for(i in idx_curr){
    g_out[[j]] <- induced_subgraph(G, V(G)$name[V(G)$name %in% names(ranked_list)[1:i]])
    j <- j+1
  }

  stat_out$n_e <- unlist(lapply(g_out, function(x) length(E(x))))
  stat_out$d <- unlist(lapply(g_out, graph.density))

  cl <- lapply(g_out, clusters)
  stat_out$n_comp <- unlist(lapply(cl, function(x) x$no))
  stat_out$n_v_max <- unlist(lapply(cl, function(x) max(x$csize)))
  stat_out$n_v_min <- unlist(lapply(cl, function(x) min(x$csize)))
  stat_out$n_v_conn_tot <- unlist(lapply(cl, function(x) sum(x$csize[x$csize>1])))
  stat_out$omega <- unlist(lapply(idx_curr, function(x) omega(G, ranked_list[1:x] / max(ranked_list[1:x]))[x]))

  return(list(sn_g=g_out, sn_summary=stat_out))


}
