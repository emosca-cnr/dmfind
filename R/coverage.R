#' Assessment of the interactome coverage of the input scores
#'
#' @param g interactome, igraph object
#' @param omic input score data on which we calculate the coverage,
#' need to have a column names "scores" (typically a p-value) and a
#' column named "gene_id" to match with the graph attribute "name"
#' @param outFile name of the output file
#' @description Assessment of the interactome coverage of the input scores
#' @usage coverage(g, omics, outFile="coverage.jpg")
#' @examples
#' \dontrun{coverage(g, omics, outFile="coverage.jpg")}
#' @import igraph
#' @export
#'

coverage <- function(g, omic, outFile="coverage.jpg") {

  omic$network <- 0
  omic$network[omic$gene_id %in% V(g)$name] <- 1
  # genes in network
  table(omic$network)

  omic <- omic[order(-omic$scores), ]
  omic$network_cov <- sapply(1:nrow(omic),
                             function(x) sum(omic$network[1:x])/x)

  png("coverage.png", width=180, height=180, units="mm", res=300)
  plot(omic$network_cov[1:5000], xlab = "score ranking", ylab= "Coverage",
       main= "Interactome coverage", type = "l", ylim=c(0,1),
       lty = 1, lwd=3, col="#158FAD")


  abline(h=mean(omic$network_cov[1:5000]), lty=2, col="#158FAD")

  dev.off()

}

