#' assess enrichment of the top networks
#' 
#' @param G igraph object
#' @param topList ranked list of vertex names
#' @param ranks ranks of topList that will be assessed
#' @param rankedVectorX0 ranked vector by X0
#' @param do.plot whether to plot or not
#' @param critical a possible rank to draw in the plot
#' @param file output file name
#' @description Perform enrichment analysis of the top network based on a ranked vector 
#' @usage assess_enrichment(G, topList, ranks, rankedVectorX0, doPlot=FALSE, critical, file="outfile.jpg")
#' @example assess_enrichment(G, topList, ranks, rankedVectorX0, doPlot=FALSE, critical, file="outfile.jpg")
#' @value randomvalue
#' @export

assess_enrichment <- function(G=NULL, topList=NULL, ranks=NULL, 
                              rankedVectorX0=NULL, doPlot=FALSE, 
                              critical=NULL, file="topNet_enrichment.jpg"){

    ranks <- seq(10, length(topList), by=10)
    gseaRes <- vector("list", length(ranks))
    netCcf <- gseaRes
    for(i in 1:length(ranks)){
        netCcf[[i]] <- get_nconn_comp(induced_subgraph(G, 
                        V(G)$name[V(G)$name %in% topList[1:ranks[i]]]), 2)

        ## if the sum of gene set scores is 0 we can not evaluate the enrichment
        if(sum(rankedVectorX0[names(rankedVectorX0) %in% 
                                V(netCcf[[i]])$name])>0){
            gseaRes[[i]] <- gsea(matrix(rankedVectorX0, ncol = 1, 
                                         dimnames = list(names(rankedVectorX0), "X0")), 
                                         list(topNet=V(netCcf[[i]])$name))
        }else{
            warning("No positive elments for ", ranks[i], "\n")
            gseaRes[[i]] <- list(gsTable=list(X0=data.frame(id=1, 
                                  es=0, pVal=1, adjpVal=1, nes=0, FDRq=1, 
                                  stringsAsFactors=FALSE)))
        }
    }

    gseaResDf <- do.call(rbind, lapply(gseaRes, function(x) x$gsTable$X0))
    gseaResDf$id <- ranks
    gseaResDf$size <- unlist(lapply(netCcf, function(x) vcount(x)))

    if(do.plot){
        jpeg(file, width=100, height=200, res=300, units="mm")
        par(mfrow=c(2, 1))
        par(mgp=c(1.5, 0.5, 0))
        par(mar=c(3, 3, 3, 1))

        top10idx <- order(log10(gseaResDf$FDRq))[1:10]

        plot(gseaResDf$id, gseaResDf$size, pch=16, xlab="rank", 
             ylab="n", cex=0.6, main = "connected components")
        if(!is.null(critical)){
            abline(v=critical, col="purple", lty=2, lwd=1)
        }
        lines(gseaResDf$id, gseaResDf$size)

        plot(gseaResDf$id, log10(gseaResDf$FDRq), pch=16, xlab="rank", 
             ylab="log10(q)", cex=0.6, main = "GSEA FDR")
        lines(gseaResDf$id, log10(gseaResDf$FDRq))
        if(!is.null(critical)){
           abline(v=critical, col="purple", lty=2, lwd=1)
        }
        abline(h=log10(0.1), col="purple", lty=2, lwd=1)
        text(gseaResDf$id[top10idx], log10(gseaResDf$FDRq)[top10idx], 
             gseaResDf$id[top10idx], cex=0.7, pos = 3)

        dev.off()
    }

    return(list(enSummary=gseaResDf, netCcf=netCcf))
}
