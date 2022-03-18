#' assess enrichment of the top networks
#' @param G igraph object
#' @param top_list ranked list of vertex names
#' @param ranks ranks of top_list that will be assessed
#' @param ranked_vector_X0 ranked vector by X0
#' @param do.plot whether to plot or not
#' @param critical a possible rank to draw in the plot
#' @param show_top_sig if not NULL, an integer that indicatesd the number of top significant ranks to display
#' @param file output file name
#'
#' @export

assess_enrichment <- function(G=NULL, top_list=NULL, ranks=NULL, ranked_vector_X0=NULL, do.plot=FALSE, critical=NULL, show_top_sig=10, file="top_net_enrichment.jpg"){

  ranks <- seq(10, length(top_list), by=10)
  gsea_res <- vector("list", length(ranks))
  net_ccf <- gsea_res
  for(i in 1:length(ranks)){
    net_ccf[[i]] <- get.n.conn.comp(induced_subgraph(G, V(G)$name[V(G)$name %in% top_list[1:ranks[i]]]), 2)

    #if the sum of gene set scores is 0 we can not evaluate the enrichment
    if(sum(ranked_vector_X0[names(ranked_vector_X0) %in% V(net_ccf[[i]])$name])>0){
      gsea_res[[i]] <- gsea(matrix(ranked_vector_X0, ncol = 1, dimnames = list(names(ranked_vector_X0), "X0")), list(top_net=V(net_ccf[[i]])$name))
    }else{
      warning("No positive elments for ", ranks[i], "\n")
      gsea_res[[i]] <- list(gs_table=list(X0=data.frame(id=1, es=0, p_val=1, adj_p_val=1, nes=0, FDRq=1, stringsAsFactors=FALSE)))
    }
  }

  gsea_res_df <- do.call(rbind, lapply(gsea_res, function(x) x$gs_table$X0))
  gsea_res_df$id <- ranks
  gsea_res_df$size <- unlist(lapply(net_ccf, function(x) vcount(x)))

  if(do.plot){
    jpeg(file, width = 180, height = 90, res=300, units="mm")
    par(mfrow=c(1, 2))
    par(mgp=c(1.5, 0.5, 0))
    par(mar=c(3, 3, 3, 1))

    if(!is.null(show_top_sig)){
      top_10_idx <- order(log10(gsea_res_df$FDRq))[1:10]
    }

    plot(gsea_res_df$id, gsea_res_df$size, pch=16, xlab="rank", ylab="n", cex=0.6, main = "connected components")
    if(!is.null(critical)){
      abline(v=critical, col="purple", lty=2, lwd=1)
    }
    lines(gsea_res_df$id, gsea_res_df$size)

    plot(gsea_res_df$id, log10(gsea_res_df$FDRq), pch=16, xlab="rank", ylab="log10(q)", cex=0.6, main = "GSEA FDR")
    lines(gsea_res_df$id, log10(gsea_res_df$FDRq))
    if(!is.null(critical)){
      abline(v=critical, col="purple", lty=2, lwd=1)
    }
    abline(h=log10(0.1), col="purple", lty=2, lwd=1)
    if(!is.null(show_top_sig)){
      text(gsea_res_df$id[top_10_idx], log10(gsea_res_df$FDRq)[top_10_idx], gsea_res_df$id[top_10_idx], cex=0.7, pos = 3)
    }


    dev.off()
  }

  return(list(en_summary=gsea_res_df, net_ccf=net_ccf))
}
