#' Plot the results of eval_eps
#'
#' @description Plot the results of eval_eps
#' @param evalEpsRes Epsilon evaluation
#' @param column which columm of evalEpsRes$csX0 and evalEpsRes$csXs
#' @param file filename
#' @examples
#' \dontrun{plot_omega_eps(evalEpsRes, file="eval_eps_X0_omega.jpg")}
#' @usage plot_omega_eps(evalEpsRes, file="eval_eps_X0_omega.jpg")
#' @return plot
#' @import grDevices
#' @import graphics
#' @import utils
#' @export

plot_omega_eps <- function(evalEpsRes=NULL, column=1, file=NULL){

    if(!is.null(file)){
        jpeg(file, width = 180, height = 100, units="mm", res=300)
    }
    par(mfrow=c(1, 2))
    par(mar=c(4, 4, 3, 3))

    plot(evalEpsRes$csX0[[1]][, column], type='l', xlab="rank", ylab="omega(X0)",
         ylim=c(0, max(unlist(lapply(evalEpsRes$csX0, function(x) x[, column])))),
         col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
    for(i in 2:length(evalEpsRes$csX0)){
        lines(evalEpsRes$csX0[[i]][, column], col=i, lty=2)
    }
    legend("bottomright", legend = paste("eps = ", format(evalEpsRes$eps, digits = 2),
           "; AUC = ", format(evalEpsRes$aucX0[, column], digits = 2), sep=""), lty=1,
           col=1:length(evalEpsRes$csX0), xpd=NA, cex=0.4)

    plot(evalEpsRes$csXs[[1]][, column], type='l', xlab="rank", ylab="omega(Xs)",
         ylim=c(0, max(unlist(lapply(evalEpsRes$csXs, function(x) x[, column])))),
         col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
    for(i in 2:length(evalEpsRes$csXs)){
        lines(evalEpsRes$csXs[[i]][, column], col=i, lty=2)
    }

    legend("bottomright", legend = paste("eps = ", format(evalEpsRes$eps, digits = 2),
          "; AUC = ", format(evalEpsRes$aucXs[, column], digits = 2), sep=""), lty=1,
          col=1:length(evalEpsRes$csXs), xpd=NA, cex=0.4)

    if(!is.null(file)){
        dev.off()
    }

}
