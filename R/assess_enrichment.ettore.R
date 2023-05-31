#' assess enrichment of the top networks
#' @param G igraph object
#' @param top_list ranked list of vertex names that will be used to define top networks
#' @param ranks ranks of top_list that will be assessed
#' @param X0_vector named numeric vector that will be tested with GSEA or ORA. In case of GSEA it will be ranked by decreasing orderg, while in the case of ORA the names of the X0 values grater than 0 will be tested for enrichment, while all X0 names will be the universe.
#' @param do.plot whether to plot or not
#' @param critical a possible rank to draw in the plot
#' @param show_top_sig if not NULL, an integer that indicatesd the number of top significant ranks to display
#' @param file output file name
#' @param k number of permutations
#'
#' @export

assess_enrichment <- function(G=NULL, top_list=NULL, ranks=NULL, X0_vector=NULL, type=c("gsea", "ora"), k=99, do.plot=FALSE, critical=NULL, show_top_sig=10, file="top_net_enrichment.jpg", min_component_size=2, min_net_size=10, min.k.nes=10){
  
  
  type <- match.arg(type, c("gsea", "ora"))
  
  if(is.null(ranks)){
    ranks <- seq(10, length(top_list), by=10)
  }
  en_res <- setNames(vector("list", length(ranks)), ranks)
  net_ccf <- en_res
  
  #defines the increasing top networks to be tested
  for(i in 1:length(ranks)){
    net_ccf[[i]] <- get.n.conn.comp(induced_subgraph(G, V(G)$name[V(G)$name %in% top_list[1:ranks[i]]]), n=min_component_size)
  }
  
  gsl <- lapply(net_ccf, function(i_net) V(i_net)$name)
  gsl <- gsl[lengths(gsl)>=min_net_size]
  net_ccf <- net_ccf[names(net_ccf) %in% names(gsl)]
  
  if(length(gsl)>0){
    
    #check X0 to choose ORA or GSEA
    if(type=="ora"){
      cat("ORA\n")
      plot_column <- "p_adj"
      
      wb <- names(X0_vector)[X0_vector>0]
      if(length(wb)>0){
        bb <- names(X0_vector)[!names(X0_vector) %in% wb]
        en_res <- ora(wb=wb, bb=bb, gsl=gsl, p_adj_method='fdr')
      }else{
        stop("Can' find any element of X0_vector > 0\n.")
      }
      
    }else{
      
      cat("GSEA\n")
      
      plot_column <- "FDRq"
      
      en_res <- gsea(rl=matrix(X0_vector, ncol = 1, dimnames = list(names(X0_vector), "X0")), gsl=gsl, k=k, min.k.nes = min.k.nes)
      en_res <- en_res$gs_table$X0
      
    }
    
    #en_res_df <- do.call(rbind, lapply(en_res, function(x) x$gs_table$X0))
    #en_res_df$id <- ranks
    en_res$size <- unlist(lapply(net_ccf, function(x) vcount(x)))
    
    if(do.plot){
      
      jpeg(file, width = 180, height = 90, res=300, units="mm")
      
      par(mfrow=c(1, 2))
      par(mgp=c(1.5, 0.5, 0))
      par(mar=c(3, 3, 3, 1))
      
      if(!is.null(show_top_sig)){
        top_10_idx <- order(log10(en_res[, plot_column]))[1:10]
      }
      
      plot(en_res$id, en_res$size, pch=16, xlab="rank", ylab="n", cex=0.6, main = "connected components")
      if(!is.null(critical)){
        abline(v=critical, col="purple", lty=2, lwd=1)
      }
      lines(en_res$id, en_res$size)
      
      plot(en_res$id, log10(en_res[, plot_column]), pch=16, xlab="rank", ylab="log10(q)", cex=0.6, main = "")
      lines(en_res$id, log10(en_res[, plot_column]))
      if(!is.null(critical)){
        abline(v=critical, col="purple", lty=2, lwd=1)
      }
      abline(h=log10(0.1), col="purple", lty=2, lwd=1)
      if(!is.null(show_top_sig)){
        text(en_res$id[top_10_idx], log10(en_res[top_10_idx, plot_column]), en_res$id[top_10_idx], cex=0.7, pos = 3)
      }
      
      
      dev.off()
    }
    
    return(list(en_summary=en_res, net_ccf=net_ccf))
    
  }else{
    return(NULL)
  }
}
