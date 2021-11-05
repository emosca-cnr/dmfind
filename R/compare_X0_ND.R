#' compare_X0_ND
#'
#' @description Generate N permutations of a binary vector
#' @param G graph
#' @param X0_ranked_by_ND ranked vector of ND values ordered by X0
#' @param X0_ranked_by_X0 ranked vector of X0 values ordered by X0
#' @param norm FALSE
#' @param do.plot FALSE
#' @return Comparison plot between X0 and ND
#' @export
#'

compare_X0_ND <- function(G=NULL, X0_ranked_by_ND=NULL, X0_ranked_by_X0=NULL, 
                          norm=FALSE, do.plot=FALSE, file="./compare_X0_ND.jpg", ...){

    ##normalization
    X0_ranked_by_ND <- X0_ranked_by_ND / max(X0_ranked_by_ND)
    X0_ranked_by_X0 <- X0_ranked_by_X0 / max(X0_ranked_by_X0)

    ans <- data.frame(
        omega_X0_X0=omega(G, u=X0_ranked_by_X0, norm=norm),
        omega_X0=omega(G, u=X0_ranked_by_ND, norm=norm),
        CCF_X0=calculate_ccf(G, names(X0_ranked_by_X0)),
        CCF=calculate_ccf(G, names(X0_ranked_by_ND)),
        row.names = 1:length(X0_ranked_by_ND))

    if(do.plot){

        jpeg(file, width=100, height=100, res=300, units="mm")
        layout(matrix(c(1, 3, 2, 3), nrow=2, byrow=T), widths=c(0.8, 0.2))
        par(mar=c(2.5, 3, 1, 1))
        par(mgp=c(1.5, 0.5, 0))

        plot(ans$omega_X0_X0, 
             ylim=c(0, max(ans[, c("omega_X0_X0", "omega_X0")])), 
             type="l", lty=2, lwd=2, xlab="rank", ylab="omega", main="")
        lines(ans$omega_X0, lwd=2)

        plot(ans$CCF_X0, ylim=c(0, max(ans[, c("CCF_X0", "CCF")])), 
             type="l", lty=2, lwd=2, xlab="rank", ylab="CCF", main="")
        lines(ans$CCF, lwd=2)

        plot.new()
        legend("right", lty=c(1, 2), legend=c("ND", "X0"), cex=0.6, xpd=T)

        dev.off()

    }

    return(ans)

}
