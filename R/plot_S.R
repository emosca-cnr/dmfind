#' Plot S
#' @description Plot S
#' @param evalEpsRes Epsilon evaluation
#' @param file filename
#' @export

plot_S <- function(calcAdjNDres=NULL, X0, initSigGenes=NULL, file="S.jpg"){

    if(!is.null(initSigGenes)){
        idxSigGenes <- which(rownames(X0) %in% initSigGenes)
    }


    jpeg(file, width = 180, height = 90, units="mm", res=300)

    par(mfrow=c(1, 2))
    par(mar=c(4, 4, 1, 1))
    plot(calcAdjNDres$S, -log10(calcAdjNDres$p), xlab="S", 
         ylab="-log10(p)", pch=16, cex=0.7)
    if(!is.null(initSigGenes)){
        points(calcAdjNDres$S[idxSigGenes], 
               -log10(calcAdjNDres$p[idxSigGenes]), col="purple", 
               cex=0.7, pch=16)
    }

    plot(X0, calcAdjNDres$Sp, xlab="X0", ylab="Sp", pch=16, cex=0.7)
    if(!is.null(initSigGenes)){
        points(X0[idxSigGenes], calcAdjNDres$Sp[idxSigGenes], 
               col="purple", cex=0.7, pch=16)
    }

    dev.off()

}
