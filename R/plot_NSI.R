#' Scatter Plot of NSI values
#' @description Plot S
#' @param nsiRes result of calc_dS() or cal_adjND()
#' @param initSigGenes NULL
#' @param file filename
#' @return plot
#' @examples
#' \dontrun{plot_S(nsiRes=NULL, X0, initSigGenes=NULL, file="S.jpg")}
#' @usage plot_S(nsiRes=NULL, X0, initSigGenes=NULL, file="S.jpg")
#' @import grDevices
#' @import graphics
#' @import utils
#' @export

plot_NSI <-
    function(nsiRes = NULL,
             x = "X0",
             y = NULL,
             column = 1,
             initSigGenes = NULL,
             file = NULL) {
        if (!is.null(initSigGenes)) {
            idxSigGenes <- which(rownames(nsiRes$X0) %in% initSigGenes)
        }
        
        
        x_data <- nsiRes[[x]][, min(column, ncol(nsiRes[[x]]))]
        y_data <- nsiRes[[y]][, min(column, ncol(nsiRes[[y]]))]
        
        if (y == "p") {
            y_data <- -log10(nsiRes[[y]][, column])
            y <- "-log10(p)"
        }
        
        if (!is.null(file)) {
            jpeg(
                file,
                width = 90,
                height = 90,
                units = "mm",
                res = 300
            )
        }
        
        par(mar = c(4, 4, 1, 1))
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
        
        if (!is.null(file)) {
            dev.off()
        }
        
    }