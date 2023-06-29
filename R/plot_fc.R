#' Plot functional cartography
#' @export

plot_fc <- function(topNetworks = NULL,
                    labelBy = NULL,
                    outFile = NULL,
                    hCut = NULL,
                    useP = FALSE) {
 
   
  wmZ <- V(topNetworks)$wmd_score
  pc <- V(topNetworks)$P
  maxZscore <- max(wmZ)
  minZscore <- min(wmZ - 1)
  
  if (is.null(hCut)) {
    if (useP) {
      hCut <- -log10(0.05)
    } else {
      hCut <- (maxZscore * 2.5) / 8
    }
  }
  
  
  if (!is.null(outFile)) {
    jpeg(
      outFile,
      width = 200,
      height = 200,
      units = "mm",
      res = 300
    )
  }
  
  layout(matrix(1:2, ncol = 2), widths = c(0.9, 0.1))
  
  par(mar = (c(4, 4, 1, 1)))
  
  plot(pc,
       wmZ,
       xlim = c(0, 1),
       xlab = "P",
       ylab = "wmd score",)
  
  #dev.off()
  ### RECT = XLEFT, YBOTTOM, XRIGHT, YTOP
  ####################################
  # horizontal cut : max z-score = author cut : author max z-score
  ### NON-HUBS: z < 2.5 - proportional to  2.5*2.5/8 = 0.78125
  
  par(xpd = F)
  
  # for R1: ULTRA PERIPHERAL NODES (p<=0.05) = nodes with all
  # their links within their module
  rect(-0.04,
       minZscore,
       0.05,
       hCut,
       col = "ivory4",
       border = NA)
  # for R2: PERIPHERAL NODES (0.05 < p <= 0.62) = nodes with
  # most links within their module
  rect(0.05,
       minZscore,
       0.62,
       hCut,
       col = "lightsalmon",
       border = NA)
  # for R3: NON-HUB CONNECTOR NODES (0.62 < p <= 0.8) = nodes
  # with many links to other modules
  rect(0.62,
       minZscore,
       0.8,
       hCut,
       col = "darkseagreen",
       border = NA)
  # for R4: NON-HUB KINLESS NODES (p > 0.8) = nodes with links
  # homogeneously distributed among all modules
  rect(0.8, minZscore, 1.3, hCut, col = "slateblue2", border = NA)
  ### HUBS : z > 2.5
  # for R5: PROVINCIAL HUBS (p <= 0.3) = hub nodes with the vast
  # majority of links within their module
  rect(-0.3,
       hCut,
       0.3,
       maxZscore + 0.2,
       col = "lightyellow",
       border = NA)
  # for R6: CONNECTOR HUBS (0.3 < p <= 0.75) = hubs with many
  # links to most of the other modules
  rect(0.3,
       hCut,
       0.75,
       maxZscore + 0.2,
       col = "mistyrose",
       border = NA)
  # for R7: KINLESS HUBS (p > 0.75) = hubs with links homogeneously
  # distributed among all modules
  rect(0.75,
       hCut,
       1.3,
       maxZscore + 0.2,
       col = "gray90",
       border = NA)
  
  points(
    pc,
    wmZ,
    xlim = c(0, 1),
    xlab = "P",
    ylab = "wmd score",
    pch = 16
  )
  
  if (!is.null(labelBy)) {
    vlab <- get.vertex.attribute(topNetworks, name = labelBy)
    thigmophobe.labels(pc, wmZ, vlab, cex = 0.8, font = 2)
  }
  
  par(mar = (c(.1, .1, 0.1, 0.1)))
  plot.new()
  legend(
    "center",
    paste0("R", 1:7),
    col = c(
      "ivory4",
      "lightsalmon",
      "darkseagreen",
      "slateblue2",
      "lightyellow",
      "mistyrose",
      "gray90"
    ),
    pch = 15,
    cex = 0.7
  )
  
  if (!is.null(outFile)) {
    dev.off()
  }
  
}