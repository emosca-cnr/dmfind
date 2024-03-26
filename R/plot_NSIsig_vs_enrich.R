#' enrichment_significance Define significance for enriched genes
#'
#' @description Define significance for enriched genes
#' @param netEnrRes dataframe resulting from 'assess_enrichment'
#' @param resSp dataframe resulting from resSp
#' @param start start
#' @param end end
#' @param by by
#' @param threshold threshold
#' @return dataframe with deltaP, deltaZ, changeRegion, euclideanDistance,
#' pc and wm_z from first and second network
#' @examples
#' \dontrun{cmp_pc_wmz(PcWmz_A, PcWmz_B)}
#' @usage cmp_pc_wmz(PcWmz_A, PcWmz_B)
#' @export
#' @importFrom plotrix thigmophobe.labels
#'

plot_NSIsig_vs_enrich <- function(netEnrRes=NULL, resSp=NULL, column=1, sigThr=0.05, sigStat = "p_adj", file=NULL) {

  p_sorted <- resSp$p[order(-resSp$Sp[, column]), column]
  
  xx <- -log10(netEnrRes$en_summary[, sigStat])
  yy <- sapply(as.numeric(netEnrRes$en_summary$rank), function(x) sum(p_sorted[1:x] < sigThr) / x)

  xx_norm <- 1 + (xx - min(xx)) / ( max(xx) - min(xx))
  yy_norm <- 1 + (yy - min(yy)) / ( max(yy) - min(yy))
  
  xyarea <- xx_norm*yy_norm
  
  if (!is.null(file)) {
    jpeg(
      file,
      width = 90,
      height = 90,
      res = 300,
      units = "mm"
    )
  }
  

  par(mar=c(3, 3, 1, 1))
  par(mgp=c(1.5, 0.5, 0))

  idx_max <- which.max(xyarea)
  
  plot(xx, yy, xlab= paste0("-log10(", sigStat, ")"), ylab=paste0("#{p < ", sigThr, "} / n"), pch=16, cex.axis=0.7, cex.lab=0.7)
  points(xx[idx_max], yy[idx_max], pch=16, col="red")
  
  text_col <- rep("black", length(xx))
  text_col[idx_max] <- "red"
  thigmophobe.labels(xx, yy, netEnrRes$en_summary$rank, cex = 0.7, col=text_col)
  
  if (!is.null(file)) {
    dev.off()
  }
  
  ans <- data.frame(id=netEnrRes$en_summary$rank, mlog10p=xx, fp=yy, mlog10p_norm=xx_norm, fp_norm=yy_norm, s=xyarea)
  
  return(ans)

}

