#' Finding significant connected components
#' @param NR_summary data.frame with NR_summary
#' @param NR_p_thr threshold for NR p value
#' @param min_rank minimum rank
#' @param max_rank maximum rank
#' @export

find_sign_conn_comp <- function(NR_summary, NR_p_thr=0.05, min_rank=50, max_rank=300){

  NR_summary_der <- NR_summary

  #NR_summary_der$rank <- 1:nrow(NR_summary)
  NR_summary_der$log10_p <- log10(NR_summary$p)
  NR_summary_der$p_der <- c(NR_summary_der$log10_p[-1] - NR_summary_der$log10_p[-length(NR_summary_der$log10_p)], 0)
  NR_summary_der$critical <- 0
  NR_summary_der$selected <- 0

  NR_summary_der$critical[NR_summary_der$p < NR_p_thr & NR_summary_der$rank >= min_rank & NR_summary_der$rank <= max_rank] <- 1

  #the critical point with the highest derivative
  critical <- NULL
  if(any(NR_summary_der$critical==1)){
    selected_rank <- NR_summary_der[NR_summary_der$critical==1, ]
    selected_rank <- selected_rank$rank[ which.max(selected_rank$p_der)]
    NR_summary_der$selected[1:selected_rank] <- 1

    critical <- max(NR_summary_der$rank[NR_summary_der$selected==1])
  }

  return(list(sign_comp_table=NR_summary_der, critical=critical))

}
