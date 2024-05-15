#' Plot the results of eval_eps
#'
#' @description Plot the results of eval_eps
#' @param evalEpsRes Epsilon evaluation
#' @param column which columm of evalEpsRes$csX0 and evalEpsRes$csXs will be used
#' @param out_dir output directory
#' @importFrom graphics legend lines par layout plot.new text barplot axis
#' @importFrom pals brewer.paired alphabet
#' @importFrom grDevices png dev.off
#' @export

plot_omega_eps <- function(evalEpsRes=NULL, column=1, out_dir=NULL){

    
    if(is.null(out_dir)){
        out_dir <- getwd()
    }else{
        dir.create(out_dir, recursive = T, showWarnings = F)
    }
    
    col <- brewer.paired(nrow(evalEpsRes$aucX0))
    
    png(file.path(out_dir, paste0("omega_rank", column, ".png")), width = 180, height = 90, res=300, units="mm")
    layout(matrix(1:3, ncol=3), widths = c(.45, .45, .1))

    par(mar=c(3, 3, .1,.1))
    par(mgp=c(2, .5, 0))
    
    plot(evalEpsRes$csX0[[1]][, column], type='l', xlab="rank", ylab="omega(X0)",
         ylim=c(0, max(unlist(lapply(evalEpsRes$csX0, function(x) x[, column])))),
         col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
    for(i in 2:length(evalEpsRes$csX0)){
        lines(evalEpsRes$csX0[[i]][, column], col=col[i], lty=2)
    }

    plot(evalEpsRes$csXs[[1]][, column], type='l', xlab="rank", ylab="omega(Xs)",
         ylim=c(0, max(unlist(lapply(evalEpsRes$csXs, function(x) x[, column])))),
         col=1, lty=2, cex.lab=0.8, cex.axis=0.8)
    for(i in 2:length(evalEpsRes$csXs)){
        lines(evalEpsRes$csXs[[i]][, column], col=col[i], lty=2)
    }

    par(mar=c(.1,.1,.1,.1))
    plot.new()
    legend("center", legend = format(evalEpsRes$eps, digits = 2), lty=1, col=col, xpd=NA, cex=0.4)
    dev.off()
    
    
    png(file.path(out_dir, paste0("omega", column, ".png")), width = 90, height = 90, res=300, units="mm")

    par(mar=c(3, 3, .1,.1))
    par(mgp=c(2, .5, 0))
    
    plot(evalEpsRes$aucX0[, column], evalEpsRes$aucXs[, column], pch=16, xlab="SUM(Omega(X0))", ylab="SUM(Omega(Xs))", cex=0.6, cex.axis=0.6, cex.lab=0.6)
    
    idx <- which(evalEpsRes$eps == evalEpsRes$opt_eps[column])
    points(evalEpsRes$aucX0[idx, column], evalEpsRes$aucXs[idx, column], pch=16, xlab="SUM(Omega(X0))", ylab="SUM(Omega(Xs))", cex=2, col="red")
    
    text(evalEpsRes$aucX0[, column], evalEpsRes$aucXs[, column], format(evalEpsRes$eps, digits = 2), cex=0.3)
    
    dev.off()
    
    png(file.path(out_dir, paste0("omega", column, ".png")), width = 90, height = 90, res=300, units="mm")
    
    par(mar=c(3, 2, 1,.1))
    par(mgp=c(1.5, .5, 0))
    
    barplot(t(cbind(X0=evalEpsRes$aucX0[, column], Xs=evalEpsRes$aucXs[, column])), names=format(evalEpsRes$eps, digits=2), xlab="eps", ylab="O(X0) + O(Xs)", cex=0.4, cex.axis=0.4, cex.lab=0.4, las=2, ylim=c(0, 2), axes=F, legend.text=T, args.legend=list(x="right", cex=0.4, bty="n"), col=alphabet(2), space=0)
    axis(2, seq(0, 2, by=.2), cex.axis=0.4, cex.lab=0.4)
    
    dev.off()
    
 
}
