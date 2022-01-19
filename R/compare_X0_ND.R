#' compare_X0_ND
#'
#' @description Generate N permutations of a binary vector
#' @param G graph
#' @param X0rankedByND ranked vector of ND values ordered by X0
#' @param X0rankedByX0 ranked vector of X0 values ordered by X0
#' @param norm FALSE
#' @param doPlot FALSE
#' @param file outfile
#' @return Comparison plot between X0 and ND
#' @examples 
#' \dontrun{compare_X0_ND(G=NULL, X0rankedByND=NULL, X0rankedByX0=NULL, norm=FALSE, doPlot=FALSE,
#'  file="compare_X0_ND.jpg")}
#' @usage compare_X0_ND(G=NULL, X0rankedByND=NULL, X0rankedByX0=NULL, norm=FALSE, doPlot=FALSE, 
#' file="compare_X0_ND.jpg")
#' @import igraph
#' @import grDevices
#' @import graphics
#' @import utils
#' @export
#'

compare_X0_ND <- function(G=NULL, X0rankedByND=NULL, X0rankedByX0=NULL, 
                          norm=FALSE, doPlot=FALSE, 
                          file="compare_X0_ND.jpg"){

    ##normalization
    X0rankedByND <- X0rankedByND / max(X0rankedByND)
    X0rankedByX0 <- X0rankedByX0 / max(X0rankedByX0)

    ans <- data.frame(
        omegaX0X0=omega(G, u=X0rankedByX0, norm=norm),
        omegaX0=omega(G, u=X0rankedByND, norm=norm),
        CCFX0=calculate_ccf(G, names(X0rankedByX0)),
        CCF=calculate_ccf(G, names(X0rankedByND)),
        row.names = 1:length(X0rankedByND))

    if(doPlot){

        jpeg(file, width=100, height=100, res=300, units="mm")
        layout(matrix(c(1, 3, 2, 3), nrow=2, byrow=TRUE), widths=c(0.8, 0.2))
        par(mar=c(2.5, 3, 1, 1))
        par(mgp=c(1.5, 0.5, 0))

        plot(ans$omegaX0X0, 
             ylim=c(0, max(ans[, c("omegaX0X0", "omegaX0")])), 
             type="l", lty=2, lwd=2, xlab="rank", ylab="omega", main="")
        lines(ans$omegaX0, lwd=2)

        plot(ans$CCFX0, ylim=c(0, max(ans[, c("CCFX0", "CCF")])), 
             type="l", lty=2, lwd=2, xlab="rank", ylab="CCF", main="")
        lines(ans$CCF, lwd=2)

        plot.new()
        legend("right", lty=c(1, 2), legend=c("ND", "X0"), cex=0.6, xpd=TRUE)

        dev.off()

    }

    return(ans)

}
