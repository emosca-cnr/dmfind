#' Plot the results of eval_eps
#'
#' @export

plot_S <- function(calc_adjND_res=NULL, X0, init_sig_genes=NULL, file="S.jpg"){


  if(!is.null(init_sig_genes)){
    idx_sig_genes <- which(rownames(X0) %in% init_sig_genes)
  }


  jpeg(file, width = 180, height = 90, units="mm", res=300)

  par(mfrow=c(1, 2))
  par(mar=c(4, 4, 1, 1))
  plot(calc_adjND_res$S, -log10(calc_adjND_res$p), xlab="S", ylab="-log10(p)", pch=16, cex=0.7)
  if(!is.null(init_sig_genes)){
    points(calc_adjND_res$S[idx_sig_genes], -log10(calc_adjND_res$p[idx_sig_genes]), col="purple", cex=0.7, pch=16)
  }

  plot(X0, calc_adjND_res$Sp, xlab="X0", ylab="Sp", pch=16, cex=0.7)
  if(!is.null(init_sig_genes)){
    points(X0[idx_sig_genes], calc_adjND_res$Sp[idx_sig_genes], col="purple", cex=0.7, pch=16)
  }

  dev.off()

}
