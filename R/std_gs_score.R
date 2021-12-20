#' Calculation of pathway score matrix
#' @param x input matrix
#' @param gsList gene set list
#' @param useMedian whether to use median (TRUE) or mean (FALSE)
#' @param scalingEps numeric value
#' @param mc.cores number of cores
#' @return matrix of scaled pathway scores
#' @examples std_gs_score(x, gsList, useMedian=FALSE, scalingEps=1, mc.cores=1)
#' @usage std_gs_score(x, gsList, useMedian=FALSE, scalingEps=1, mc.cores=1)
#' @export
#'
std_gs_score <- function(x, gsList, useMedian=FALSE, scalingEps=1, mc.cores=1){

    #calculation of gsScore and scaled gsScore
    out <- split(x, rep(1:ncol(x), each = nrow(x)))

    if(mc.cores > 1){
        out <- as.data.frame(BiocParallel::bplapply(out, gsScore, 
                        gsList=gsList, mc.cores=mc.cores))
    }else{
        out <- as.data.frame(lapply(out, gsScore, gsList=gsList))
    }

    colnames(out) <- colnames(x)

    if(useMedian){
        out <- t(scale(t(out), center=apply(out, 1, median), 
                       scale=apply(out, 1, sd) + scalingEps))
    }else{
        out <- t(scale(t(out), center=apply(out, 1, mean), 
                       scale=apply(out, 1, sd) + scalingEps))
    }

    return(out)

}
