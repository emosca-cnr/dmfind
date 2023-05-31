#' Plot the results of eval_eps
#' 
#' @description Plot the results of eval_eps
#' @param nrRes network resampling result
#' @param signCompTable ?
#' @param file filename
#' @return plot
#' @examples 
#' \dontrun{plot_NR(nrRes=NULL, signCompTable=NULL, file="NR.jpg")}
#' @usage plot_NR(nrRes=NULL, signCompTable=NULL, file="NR.jpg")
#' @import grDevices
#' @import graphics
#' @import utils
#' @export

plot_NR <- function(nrRes=NULL, signCompTable=NULL, file="NR.jpg"){

    if(!is.null(signCompTable)){
        selectedRanks <- max(signCompTable$rank[signCompTable$selected==1])
    }

    jpeg(file, width = 180, height = 90, units="mm", res=300)

    par(mfrow=c(1, 2))
    par(mar=c(4, 4, 1, 1))
    plot(nrRes$omegaPerm[, 1], xlab="rank", ylab="omega", type="l", col="gray",
         ylim = c(0, max(c(nrRes$NRsummary$omega, unlist(nrRes$omegaPerm)))), 
         lty=2, lwd=2)
    for(i in 1:min(100, ncol(nrRes$omegaPerm))){
        lines(nrRes$omegaPerm[, i], col="gray", lty=2)
    }
    lines(nrRes$NRsummary$omega, col="red")

    if(!is.null(signCompTable)){
        abline(v=selectedRanks, col="purple", lty=2, lwd=1)
    }

    plot(log10(nrRes$NRsummary$p), xlab="rank", ylab="log10(p)", type="l")
    points(log10(nrRes$NRsummary$p), pch=16, cex=0.5)
    if(!is.null(signCompTable)){
        idxCritical <- which(signCompTable$critical==1)
        points(nrRes$NRsummary$rank[idxCritical], 
               log10(nrRes$NRsummary$p[idxCritical]), pch=16, 
               cex=0.5, col="purple")
        abline(v=selectedRanks, col="purple", lty=2, lwd=1)
    }

    dev.off()

}
