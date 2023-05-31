#' enrichment_significance Define significance for enriched genes
#'
#' @description Define significance for enriched genes
#' @param ae_res_full dataframe resulting from 'assess_enrichment'
#' @param mNDp dataframe resulting from mNDp
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
#'
#'


enrichment_significance <- function(ae_res_full, mNDp, start, end, by, threshold) {

  library(plotrix)

  mNDp <- mNDp[order(-mNDp$mNDp), ]
  xx <- -log10(ae_res_full$en_summary$FDRq)
  yy <- sapply(seq(50, 500, by=10), function(x) sum(res_ann$p[1:x] < 0.005) / x)

  jpeg("enrichment_significance.jpg", width = 90, height = 90, res=300, units="mm")

  par(mar=c(3, 3, 0.1, 1))
  par(mgp=c(1.5, 0.5, 0))

  plot(xx, yy, xlab="-log10(q)", ylab="#{p < 0.05} / n", pch=16, cex.axis=0.7, cex.lab=0.7)
  points(xx[36], yy[36], pch=16, col="red")
  thigmophobe.labels(xx, yy, seq(50, 500, by=10), cex=0.6, xpd=T)

  dev.off()

  df <- data.frame(xx, yy)
  return(df)

}

