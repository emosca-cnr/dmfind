#' Pipeline for the extraction of differentially enriched modules from an input vector
#' @param X input vector
#' @param mc.cores number of cores
#' @param G gene x gene undirected interaction graph
#' @param NR_k number of permutations
#' @param NR_p_thr threshold for NR p value
#' @param min_rank minimum rank
#' @param max_rank maximum rank
#' @param min_module_size_final minimun size of final module
#' @param min_subnet_size minimum number of genes connected in the subnetworks to be extracted
#' @param plot_outfile file name of the output figure
#' @param NR_max_rank network resampling will be applied from 1 to NR_max_rank
#' @param plot_flag whether to plot or not the network
#' @param ... further arguments to igraph.plot
#' @import parallel
#' @export
#'

sv2dem <- function(X, G, mc.cores=2, NR_p_thr=0.05, NR_k=100, min_module_size_final=10, min_rank=100, max_rank=500, NR_max_rank=600, min_subnet_size=2, plot_flag=FALSE, plot_outfile="graph.jpg", ...){

  if(min_rank > length(X)){
    stop("min_rank can't be larger than X")
  }
  if(max_rank > length(X)){
    stop("max_rank can't be larger than X")
  }
  if(min_module_size_final > length(X)){
    stop("min_module_size_final can't be larger than X")
  }
  if(NR_max_rank > length(X)){
    warning("adjusting NR_max_rank to lenght(X)")
    NR_max_rank <- min(NR_max_rank, length(X))
  }

  cat('assembling smoothed scores\n')
  smoothed_scores <- data.frame(id=names(X), X=X, stringsAsFactors = FALSE)

  ranked_vector <- sort(array(smoothed_scores$X, dimnames = list(smoothed_scores$id)), decreasing = TRUE)[1:NR_max_rank]

  cat('NR\n')
  NR_res <- NR(G, ranked_vector, NR_k, mc.cores)

  cat('find_sign_conn_comp\n')
  sig_comp <- find_sign_conn_comp(NR_res$NR_summary, NR_p_thr, min_rank, max_rank)

  summary_log <- list(length=length(X))

  if(sum(sig_comp$selected) >= min_module_size_final){
    cat('extract_module\n')
    modules <- extract_module(G, sig_comp$id[sig_comp$selected==1], vertices_weight=ranked_vector, min_subnet_size = min_subnet_size, plot_outfile=plot_outfile, ...)
    summary_log$modules_size <- length(V(modules$conn_subnets))
    summary_log$modules_links <- length(E(modules$conn_subnets))
    summary_log$modules_d <- graph.density(modules$conn_subnets)
  }else{
    cat("Can't find any module of size >=", min_module_size_final, "\n")
    modules <- NA
  }

  settings <- list(mc.cores=mc.cores, NR_p_thr=NR_p_thr, NR_k=NR_k, min_module_size_final=min_module_size_final, min_rank=min_rank, max_rank=max_rank, NR_max_rank=NR_max_rank, min_subnet_size=min_subnet_size, plot_outfile=plot_outfile)

  return(list(smoothed_scores=smoothed_scores, NR_res=NR_res, sig_comp=sig_comp, modules=modules, summary_log=summary_log, settings=settings))

}
