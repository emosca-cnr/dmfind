#' Plot network enrichment results
#'
#' @param showTopSig if not NULL, an integer that indicatesd the number of top significant ranks to display
#' @param sigStat the quantity to be plotted: for ORA, can be any of "er", "p" or "p_adj"; for GSEA can be "es", "p", "p_val", "nes" or "FDRq"
#' @export

plot_net_enrich <-
  function(netEnrRes = NULL,
           showTopSig = 10,
           sigStat = "p_adj",
           file = NULL) {
    en_res <- netEnrRes$en_summary
    
    if (!is.null(file)) {
      jpeg(
        file,
        width = 180,
        height = 90,
        res = 300,
        units = "mm"
      )
    }
    
    
    par(mfrow = c(1, 2))
    par(mgp = c(1.5, 0.5, 0))
    par(mar = c(3, 3, 3, 1))
    
    if (sigStat %in% c("p", "p_adj", "FDRq")) {
      yy <- log10(en_res[, sigStat])
      ylab <- paste0("log10(", sigStat, ")")
    } else{
      yy <- en_res[, sigStat]
      ylab <- sigStat
    }
    
    top_10_idx <-
      order(yy, decreasing = !(sigStat %in% c("p", "p_adj", "FDRq")))[1:showTopSig]
    
    plot(
      en_res$rank,
      en_res$size,
      pch = 16,
      xlab = "rank",
      ylab = "n",
      cex = 0.6,
      main = "connected components"
    )
    lines(en_res$rank, en_res$size)
    
    plot(
      en_res$rank,
      yy,
      pch = 16,
      xlab = "rank",
      ylab = ylab,
      cex = 0.6,
      main = "enrichment"
    )
    lines(en_res$rank, yy)
    
    if (!is.null(showTopSig)) {
      points(en_res$rank[top_10_idx],
             yy[top_10_idx],
             pch = 16,
             cex = 0.6,
             col = "red")
      thigmophobe.labels(en_res$rank[top_10_idx], yy[top_10_idx], en_res$rank[top_10_idx], cex =
                           0.7)#, pos = 3)
    }
    
    
    if (!is.null(file)) {
      dev.off()
    }
  }