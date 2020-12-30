#' Pipeline for the extraction of differentially enriched modules from smoothed matrices
#' @param X0 input matrix
#' @param Xs input matrix
#' @param ranking_type criterion for ranking the molecular entities
#' @param classes numeric vector of \{1,2\}
#' @param eps numeric value
#' @param use.exp TRUE/FALSE whether to use the exponent at the denominator
#' @param median TRUE/FALSE whether to use median instead of mean
#' @param classes_perm list of numeric vectors of \{1,2\}
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
#' @param vector_mode TRUE/FALSE
#' @param ... further arguments to igraph.plot
#' @import parallel
#' @export
#'

sm2dem <- function(X0, Xs, G, classes, eps=1, use.exp=FALSE, median=FALSE, classes_perm=NULL, mc.cores=2, ranking_type=c("dR", "dS", "dXs", "dX0"), NR_p_thr=0.05, NR_k=100, min_module_size_final=10, min_rank=100, max_rank=500, NR_max_rank=600, min_subnet_size=2, plot_flag=FALSE, plot_outfile="graph.jpg", vector_mode=FALSE, ...){

  ranking_type <- match.arg(ranking_type)

  cat('calc_dS\n')
  calc_dS_res <- calc_dS(X0, Xs, classes, eps=eps, use.exp=use.exp, median=median, classes_perm=classes_perm, mc.cores=mc.cores, vector_mode=vector_mode)

  if(ranking_type == 'dR'){
    ranked_vector <- sort(array(calc_dS_res$dR, dimnames = list(calc_dS_res$id)), decreasing = TRUE)[1:NR_max_rank]
  }
  if(ranking_type == 'dS'){
    ranked_vector <- sort(array(calc_dS_res$dS, dimnames = list(calc_dS_res$id)), decreasing = TRUE)[1:NR_max_rank]
  }
  if(ranking_type == 'dX0'){
    ranked_vector <- sort(array(calc_dS_res$dX0, dimnames = list(calc_dS_res$id)), decreasing = TRUE)[1:NR_max_rank]
  }
  if(ranking_type == 'dXs'){
    ranked_vector <- sort(array(calc_dS_res$dXs, dimnames = list(calc_dS_res$id)), decreasing = TRUE)[1:NR_max_rank]
  }


  cat('NR\n')
  NR_res <- NR(G, ranked_vector, NR_k, mc.cores)

  cat('find_sign_conn_comp\n')
  sig_comp <- find_sign_conn_comp(NR_res$NR_summary, NR_p_thr, min_rank, max_rank)

  summary_log <- list(rows=nrow(X0), cols=ncol(X0))

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

  settings <- list(classes, eps=eps, use.exp=use.exp, median=median, mc.cores=mc.cores, ranking_type=ranking_type, NR_p_thr=NR_p_thr, NR_k=NR_k, min_module_size_final=min_module_size_final, min_rank=min_rank, max_rank=max_rank, NR_max_rank=NR_max_rank, min_subnet_size=min_subnet_size, plot_outfile=plot_outfile)

  return(list(calc_dS_res=calc_dS_res, NR_res=NR_res, sig_comp=sig_comp, modules=modules, summary_log=summary_log, settings=settings))

}
