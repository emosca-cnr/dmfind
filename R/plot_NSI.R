#' Scatter Plot of NSI values
#' @description Plot NSI
#' @param nsiRes result of calc_dS() or cal_adjND()
#' @param x what to put on the horizontal axis
#' @param y what to put on the vertical axis
#' @param initSigGenes optional vector of initial sifnigicant genes
#' @param column column to consider
#' @importFrom graphics par points legend
#' @export

plot_NSI <-
    function(nsiRes = NULL, x = "X0", y = NULL, column = 1, initSigGenes = NULL) {
        
        if (!is.null(initSigGenes)) {
            idxSigGenes <- which(rownames(nsiRes$X0) %in% initSigGenes)
        }
        
        x_data <- nsiRes[[x]][, min(column, ncol(nsiRes[[x]]))]
        y_data <- nsiRes[[y]][, min(column, ncol(nsiRes[[y]]))]
        
        if (y %in% c("p", "pX0", "pW")) {
            y_data <- -log10(nsiRes[[y]][, column])
            y <- paste0("-log10(", y, ")")
        }
        
        plot(
            x_data,
            y_data,
            xlab = x,
            ylab = y,
            pch = 16,
            cex = 0.7,
            cex.axis = 0.7,
            cex.lab = 0.7
        )
        if (!is.null(initSigGenes)) {
            points(x_data[idxSigGenes],
                   y_data[idxSigGenes],
                   col = "purple",
                   cex = 0.7,
                   pch = 16)
            legend("bottomright", c("high X0", "others"), pch=16, col=c("purple", "black"), cex=0.7)
        }
        
    }
