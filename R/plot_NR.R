#' Plot the results of eval_eps
#'
#' @export

plot_NR <- function(nr_res=NULL, sign_comp_table=NULL, file="NR.jpg"){


  if(!is.null(sign_comp_table)){
    selected_ranks <- max(sign_comp_table$rank[sign_comp_table$selected==1])
  }

  jpeg(file, width = 180, height = 90, units="mm", res=300)

  par(mfrow=c(1, 2))
  par(mar=c(4, 4, 1, 1))
  plot(nr_res$omega_perm[, 1], xlab="rank", ylab="omega", type="l", col="gray", ylim = c(0, max(c(nr_res$NR_summary$omega, unlist(nr_res$omega_perm)))), lty=2, lwd=2)
  for(i in sample(2:99, 40)){
    lines(nr_res$omega_perm[, i], col="gray", lty=2)
  }
  lines(nr_res$NR_summary$omega, col="red")

  if(!is.null(sign_comp_table)){
    abline(v=selected_ranks, col="purple", lty=2, lwd=1)
  }


  plot(log10(nr_res$NR_summary$p), xlab="rank", ylab="log10(p)", type="l")
  points(log10(nr_res$NR_summary$p), pch=16, cex=0.5)
  if(!is.null(sign_comp_table)){
    idx_critical <- which(sign_comp_table$critical==1)
    points(nr_res$NR_summary$rank[idx_critical], log10(nr_res$NR_summary$p[idx_critical]), pch=16, cex=0.5, col="purple")
    abline(v=selected_ranks, col="purple", lty=2, lwd=1)
  }

  dev.off()

}
