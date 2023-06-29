#' Compare top networks scores
#'
#' @description Define significance for enriched genes
#' @param ae_res_full dataframe resulting from 'assess_enrichment'
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

cmp_top_net_scores <-
  function(NRRes = NULL,
           netEnrRes = NULL,
           column = 1,
           sigStatNR = "p",
           sigStatEn = "p",
           norm = TRUE,
           top=5,
           file = NULL) {
    #p_sorted <- resSp$p[order(-resSp$Sp[, column]), column]
    
    ans <-
      merge(
        NRRes$NRsummary[, c("rank", sigStatNR)],
        netEnrRes$en_summary[, c("rank", sigStatEn)],
        by = "rank",
        suffixes = c("_NR", "_NE"),
        sort = F
      )
    ans <- ans[order(as.numeric(ans$rank)),]
    
    xx <- -log10(ans[, 2])
    #yy <- sapply(as.numeric(netEnrRes$en_summary$id), function(x) sum(p_sorted[1:x] < sigThr) / x)
    yy <- -log10(ans[, 3])
    
    if (norm) {
      xx <- 1 + (xx - min(xx)) / (max(xx) - min(xx))
      yy <- 1 + (yy - min(yy)) / (max(yy) - min(yy))
    }
    
    #xyarea <- xx_norm*yy_norm
    ans$score <- xx * yy
    
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
    par(mar = c(3, 3, 1, 1))
    par(mgp = c(1.5, 0.5, 0))
    
    idx_max <- order(-ans$score)[1:top]
    
    if (nrow(unique(ans[, 2:3])) != nrow(ans)) {
      cat("Overlapping points detected, adding some jitter\n")
      xx <- jitter(xx)
      yy <- jitter(yy)
    }
    
    plot(
      xx,
      yy,
      xlab = paste0("NR [-log10(", sigStatNR, ")]"),
      ylab = paste0("NE [-log10(", sigStatEn, ")]"),
      pch = 16,
      cex = 0.6,
      cex.axis = 0.7,
      cex.lab = 0.7
    )
    points(xx[idx_max],
           yy[idx_max],
           pch = 16,
           col = "red",
           cex = 0.6)
    
    text_col <- rep("black", length(xx))
    text_col[idx_max] <- "red"
    thigmophobe.labels(xx, yy, ans$rank, cex = 0.6, col = text_col)
    
    plot(
      as.numeric(ans$rank),
      ans$score,
      xlab = "rank",
      ylab = "[-log10(p)]^2",
      pch = 16,
      cex = 0.6,
      cex.axis = 0.7,
      cex.lab = 0.7,
      col=text_col
    )
    
    if (!is.null(file)) {
      dev.off()
    }
    
    return(ans)
    
  }
