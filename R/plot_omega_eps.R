#' Plot the results of eval_eps
#'
#' @export

plot_omega_eps <- function(eval_eps_res, file="eval_eps_X0_omega.jpg"){

  jpeg(file, width = 180, height = 100, units="mm", res=300)
  par(mfrow=c(1, 2))
  par(mar=c(4, 4, 3, 3))

  plot(eval_eps_res$csX0[[1]][, 1], type='l', xlab="rank", ylab="omega(X0)", ylim=c(0, max(unlist(lapply(eval_eps_res$csX0, function(x) x[, 1])))), col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
  for(i in 2:nrow(eps)){
    lines(eval_eps_res$csX0[[i]][, 1], col=i, lty=2)
  }
  legend("topleft", legend = paste("eps = ", format(eps[, 1], digits = 2), "; AUC = ", format(eval_eps_res$aucX0, digits = 2), sep=""), lty=1, col=1:nrow(eps), xpd=NA, cex=0.4)

  plot(eval_eps_res$csXs[[1]][, 1], type='l', xlab="rank", ylab="omega(Xs)", ylim=c(0, max(unlist(lapply(eval_eps_res$csXs, function(x) x[, 1])))), col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
  for(i in 2:nrow(eps)){
    lines(eval_eps_res$csXs[[i]][, 1], col=i, lty=2)
  }

  legend("topleft", legend = paste("eps = ", format(eps[, 1], digits = 2), "; AUC = ", format(eval_eps_res$aucXs, digits = 2), sep=""), lty=1, col=1:nrow(eps), xpd=NA, cex=0.4)

  dev.off()

}
