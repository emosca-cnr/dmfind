#' Plot the results of eval_eps
#'
#' @description Plot the results of eval_eps
#' @param evalEpsRes Epsilon evaluation
#' @param file filename
#' @examples plot_omega_eps(evalEpsRes, file="eval_eps_X0_omega.jpg")
#' @usage plot_omega_eps(evalEpsRes, file="eval_eps_X0_omega.jpg")
#' @return plot
#' @export

plot_omega_eps <- function(evalEpsRes, file="eval_eps_X0_omega.jpg"){

    jpeg(file, width = 180, height = 100, units="mm", res=300)
    par(mfrow=c(1, 2))
    par(mar=c(4, 4, 3, 3))

    plot(evalEpsRes$csX0[[1]][, 1], type='l', xlab="rank", ylab="omega(X0)", 
         ylim=c(0, max(unlist(lapply(evalEpsRes$csX0, function(x) x[, 1])))), 
         col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
    for(i in 2:nrow(eps)){
        lines(evalEpsRes$csX0[[i]][, 1], col=i, lty=2)
    }
    legend("topleft", legend = paste("eps = ", format(eps[, 1], digits = 2), 
           "; AUC = ", format(evalEpsRes$aucX0, digits = 2), sep=""), lty=1, 
           col=1:nrow(eps), xpd=NA, cex=0.4)

    plot(evalEpsRes$csXs[[1]][, 1], type='l', xlab="rank", ylab="omega(Xs)", 
         ylim=c(0, max(unlist(lapply(evalEpsRes$csXs, function(x) x[, 1])))), 
         col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
    for(i in 2:nrow(eps)){
        lines(evalEpsRes$csXs[[i]][, 1], col=i, lty=2)
    }

    legend("topleft", legend = paste("eps = ", format(eps[, 1], digits = 2), 
          "; AUC = ", format(evalEpsRes$aucXs, digits = 2), sep=""), lty=1, 
          col=1:nrow(eps), xpd=NA, cex=0.4)

    dev.off()

}
