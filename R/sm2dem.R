#' Pipeline for the extraction of differentially enriched modules from smoothed matrices
#' @param X0 input matrix
#' @param Xs input matrix
#' @param rankingType criterion for ranking the molecular entities
#' @param classes numeric vector of \{1,2\}
#' @param eps numeric value
#' @param useExp TRUE/FALSE whether to use the exponent at the denominator
#' @param median TRUE/FALSE whether to use median instead of mean
#' @param classesPerm list of numeric vectors of \{1,2\}
#' @param mc.cores number of cores
#' @param G gene x gene undirected interaction graph
#' @param NRk number of permutations
#' @param NRpThr threshold for NR p value
#' @param minRank minimum rank
#' @param maxRank maximum rank
#' @param minModuleSizeFinal minimun size of final module
#' @param minSubnetSize minimum number of genes connected in the subnetworks to be extracted
#' @param plotOutfile file name of the output figure
#' @param NRmaxRank network resampling will be applied from 1 to NRmaxRank
#' @param plotFlag whether to plot or not the network
#' @param vectorMode TRUE/FALSE
#' @param ... further arguments to igraph.plot
#' @import BiocParallel
#' @export
#'

sm2dem <- function(X0, Xs, G, classes, eps=1, useExp=FALSE, median=FALSE, 
                   classesPerm=NULL, mc.cores=2, 
                   rankingType=c("dR", "dS", "dXs", "dX0"), 
                   NRpThr=0.05, NRk=100, minModuleSizeFinal=10, 
                   minRank=100, maxRank=500, NRmaxRank=600, minSubnetSize=2, 
                   plotFlag=FALSE, plotOutfile="graph.jpg", 
                   vectorMode=FALSE, ...){

    rankingType <- match.arg(rankingType)

    cat('calc_dS\n')
    calcdSres <- calc_dS(X0, Xs, classes, eps=eps, useExp=useExp, 
                         median=median, classesPerm=classesPerm, 
                         mc.cores=mc.cores, vectorMode=vectorMode)

    if(rankingType == 'dR'){
        rankedVector <- sort(array(calcdSres$dR, 
                             dimnames = list(calcdSres$id)), 
                             decreasing = TRUE)[1:NRmaxRank]
    }
    if(rankingType == 'dS'){
        rankedVector <- sort(array(calcdSres$dS, 
                             dimnames = list(calcdSres$id)), 
                             decreasing = TRUE)[1:NRmaxRank]
    }
    if(rankingType == 'dX0'){
        rankedVector <- sort(array(calcdSres$dX0, 
                             dimnames = list(calcdSres$id)), 
                             decreasing = TRUE)[1:NRmaxRank]
    }
    if(rankingType == 'dXs'){
        rankedVector <- sort(array(calcdSres$dXs, 
                             dimnames = list(calcdSres$id)), 
                             decreasing = TRUE)[1:NRmaxRank]
    }

    cat('NR\n')
    NRres <- NR(G, rankedVector, NRk, mc.cores)

    cat('find_sign_conn_comp\n')
    sigComp <- find_sign_conn_comp(NRres$NR_summary, NRpThr, minRank, maxRank)

    summaryLog <- list(rows=nrow(X0), cols=ncol(X0))

    if(sum(sigComp$selected) >= minModuleSizeFinal){
        cat('extract_module\n')
        modules <- extract_module(G, sigComp$id[sigComp$selected==1],
                              vertices_weight=rankedVector, 
                              minSubnetSize = minSubnetSize, 
                              plotOutfile=plotOutfile, ...)
        summaryLog$modulesSize <- length(V(modules$connSubnets))
        summaryLog$modulesLinks <- length(E(modules$connSubnets))
        summaryLog$modulesD <- graph.density(modules$connSubnets)
    }else{
        cat("Can't find any module of size >=", minModuleSizeFinal, "\n")
        modules <- NA
    }

    settings <- list(classes, eps=eps, useExp=useExp, median=median, 
                   mc.cores=mc.cores, rankingType=rankingType, NRpThr=NRpThr, 
                   NRk=NRk, minModuleSizeFinal=minModuleSizeFinal, 
                   minRank=minRank, maxRank=maxRank, NRmaxRank=NRmaxRank, 
                   minSubnetSize=minSubnetSize, plotOutfile=plotOutfile)

    return(list(calcdSres=calcdSres, NRres=NRres, sigComp=sigComp, 
              modules=modules, summaryLog=summaryLog, settings=settings))

}
