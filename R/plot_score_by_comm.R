#' plot_score_by_comm
#' @param g igraph object
#' @param comm_attr attrubute that contains community information
#' @param score_attr  attrubute that contains gene scores; high values are considered significant
#' @param min_comm_size minimum community size
#' @param top_genes_to_label number of genes to label
#' @param out_file out file name
#' @importFrom plotrix thigmophobe.labels
#' @import igraph
#' @export

plot_score_by_comm <- function(g=NULL, comm_attr="comm_id", score_attr=NULL, min_comm_size=5, top_genes_to_label=3, out_file="score_by_comm.png"){
  
  comm_id <- get.vertex.attribute(graph = g, name = comm_attr)
  score_val <- get.vertex.attribute(graph = g, name = score_attr)
  
  png(out_file, width=180, height=180, res=300, units="mm")
  
  comm_order <- comm_id
  comm_order <- comm_order[comm_order >= min_comm_size]
  
  top_genes_to_label <- 3
  
  i<-1
  xx <- jitter(rep(i, comm_order[i]), amount = 0.2)
  yy <- score_val[comm_id == names(comm_order)[i]]

  plot(xx, yy, xlim = c(0, length(comm_order)+1), ylim = c(0, max(score_val)), xlab="# community", ylab="y", pch=16, cex=0.7, xaxt="n")
  axis(1, at=c(1:length(comm_order)), labels = names(comm_order))
  boxplot(yy, at=i, add=T, outline=F, col=NA)
  
  if(any(yy>0.05)){
    top_gg <- names(yy)[order(-yy)[1:min(top_genes_to_label, comm_order[i])]]
    plotrix::thigmophobe.labels(xx[1:length(top_gg)], yy[names(yy) %in% top_gg], V(g)$label[match(top_gg, V(g)$name)], cex=0.7, font=2)
  }
  
  for(i in 2:length(comm_order)){
    
    xx <- jitter(rep(i, comm_order[i]), amount = 0.2)
    yy <- score_val[comm_id == names(comm_order)[i]]

    points(xx, yy, pch=16, cex=0.7)
    boxplot(yy, at=i, add=T, outline = F, col=NA)
    
    if(any(yy>0.05)){
      top_gg <- names(yy)[order(-yy)[1:min(top_genes_to_label, comm_order[i])]]
      plotrix::thigmophobe.labels(xx[1:length(top_gg)], yy[names(yy) %in% top_gg], V(g)$label[match(top_gg, V(g)$name)], cex=0.7, font=2)
    }
  }
  
  dev.off()
  
}
