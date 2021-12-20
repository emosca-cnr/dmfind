#' pc_wmz Calculate participation coeffcient and within module degree z-score
#' 
#' @description Calculate participation coeffcient and within module degree z-score
#' @param topNetwork network with some parameters: "comm_id" for community membership and "label" for labels (for plotting purpose).
#' @param allNames TRUE/FALSE. Used to plot all gene names.
#' @param extremePoints TRUE/FALSE. Used to plot extreme poits of the plot, such as the gene with max and min participation coefficient and z-score.
#' @param chooseYourGenes Vector of names (matching with "label" of network). Used to plot genes name of your choice.
#' @return dataframe with pc, wmz, region, description and role
#' @examples pc_wmz(topNetwork, allNames = FALSE, extremePoints = FALSE, chooseYourGenes = NULL, outFile="pc_wmz.jpg", doPlot=TRUE)
#' @usage pc_wmz(topNetwork, allNames = FALSE, extremePoints = FALSE, chooseYourGenes = NULL, outFile="pc_wmz.jpg", doPlot=TRUE)
#' @import brainGraph plotrix
#' @export

pc_wmz <- function(topNetwork, allNames=FALSE, extremePoints=FALSE, chooseYourGenes=NULL, outFile="pc_wmz.jpg", doPlot=TRUE) {

    # partecipation coefficient
    pc <- brainGraph::part_coeff(topNetwork, V(topNetwork)$comm_id)

    # within module z-score
    wmZ <-brainGraph::within_module_deg_z_score(topNetwork, V(topNetwork)$comm_id)
    wmZ[is.na(wmZ)] <- 0

    # cbind pc and wmZ (x and y coordinates)
    df <- data.frame(name=names(pc), label=V(topNetwork)$label[match(names(pc), 
                    V(topNetwork)$name)], comm_id=V(topNetwork)$comm_id[match(names(pc),
                    V(topNetwork)$name)], pc=pc, wmZ=wmZ, stringsAsFactors = FALSE)

    maxZscore <- max(wmZ)
    hCut <- (maxZscore * 2.5) / 8

    # look for coordinates
    # NON-HUBS : z < 2.5
    # for R1: ULTRA PERIPHERAL NODES (p<=0.05) = nodes with all their links within their module
    df$region[df$pc <= 0.05 & df$wmZ < hCut] <- "R1"
    df$definition[df$region == "R1"] <- "ULTRA PERIPHERAL NODES"
    df$description[df$region == "R1"] <- 
        "nodes with all their links within their module"

    # for R2: PERIPHERAL NODES (0.05 < p <= 0.62) = nodes with most links within their module
    df$region[df$pc > 0.05 & df$pc <= 0.62 & df$wmZ < hCut] <- "R2"
    df$definition[df$region == "R2"] <- "PERIPHERAL NODES"
    df$description[df$region == "R2"] <- 
        "nodes with most links within their module"

    # for R3: NON-HUB CONNECTOR NODES (0.62 < p <= 0.8) = nodes with many links to other modules
    df$region[df$pc > 0.62 & df$pc <= 0.8 & df$wmZ < hCut] <- "R3"
    df$definition[df$region == "R3"] <- "NON-HUB CONNECTOR NODES"
    df$description[df$region == "R3"] <- 
        "nodes with many links to other modules"

    # for R4: NON-HUB KINLESS NODES (p > 0.8) = nodes with links homogeneously distributed among all modules
    df$region[df$pc > 0.8 & df$wmZ < hCut] <- "R4"
    df$definition[df$region == "R4"] <- "NON-HUB KINLESS NODES"
    df$description[df$region == "R4"] <- 
        "nodes with links homogeneously distributed among all modules"

    ### HUBS : z > 2.5
    # for R5: PROVINCIAL HUBS (p <= 0.3) = hub nodes with the vast majority of links within their module
    df$region[df$pc <= 0.3 & df$wmZ > hCut] <- "R5"
    df$definition[df$region == "R5"] <- "PROVINCIAL HUBS"
    df$description[df$region == "R5"] <- 
        "hub nodes with the vast majority of links within their module"

    # for R6: CONNECTOR HUBS (0.3 < p <= 0.75) = hubs with many links to most of the other modules
    df$region[df$pc > 0.3 & df$pc <= 0.75 & df$wmZ > hCut] <- "R6"
    df$definition[df$region == "R6"] <- "CONNECTOR HUBS"
    df$description[df$region == "R6"] <- 
        "hubs with many links to most of the other modules"

    # for R7: KINLESS HUBS (p > 0.75) = hubs with links homogeneously distributed among all modules
    df$region[df$pc > 0.75 & df$wmZ > hCut] <- "R7"
    df$definition[df$region == "R7"] <- "KINLESS HUBS"
    df$description[df$region == "R7"] <- 
        "hubs with links homogeneously distributed among all modules"

    ########################################
    ########### check parameters ###########
    ########################################

    if(doPlot){
        if (allNames == TRUE & extremePoints == TRUE | allNames == TRUE & 
            extremePoints == TRUE & !is.null(chooseYourGenes)) {
        stop("Impossibile to satisfy all the conditions. 
            Choose different parameters.")
        } else {

    ########################################
    ################# PLOT #################
    ########################################

    vLabl <- ""

    if (allNames == TRUE) {
        warning("You have selected to plot all gene names. 
                To plot only expreme points use 'extremePoints = TRUE")
        vLabl <- V(topNetwork)$label
    }
    if (extremePoints == TRUE) {
        warning("You have selected to plot extreme gene names. 
                To plot all the gene names use 'allNames = TRUE'")

        nameMaxZscore <- 
            V(topNetwork)$label[match(names(pc)[which.max(wmZ)], 
                                       V(topNetwork)$name)]

        # gene with max pc
        nameMaxPc <- 
            V(topNetwork)$label[match(names(pc)[which.max(pc)], 
                                       V(topNetwork)$name)]

        # genes with pc = 0
        V(topNetwork)$label[match(names(pc)[pc == 0 & which.max(wmZ)], 
                                   V(topNetwork)$name)]

        # gene with min z score
        nameMinZscore <- 
            V(topNetwork)$label[match(names(pc)[which.min(wmZ)], 
                                       V(topNetwork)$name)]

        # gene with min pc
        nameMinPc <- 
            V(topNetwork)$label[match(names(pc)[which.min(pc)], 
                                       V(topNetwork)$name)]

        vLabl <- V(topNetwork)$label
        vLabl[!vLabl %in% c(nameMaxZscore, 
                              nameMaxPc, nameMinZscore, nameMinPc)] <- ""
      }

      if (!is.null(chooseYourGenes)) {

        warning("You have selected to plot your gene names. To plot all the gene names use 'allNames = TRUE'")
        vLabl <- V(topNetwork)$label
        vLabl[!vLabl %in% chooseYourGenes] <- ""

      }

      #V(topNetwork)$vLabl <- vLabl

      # PLOT ALL GENE NAMES
      if (allNames == FALSE & extremePoints == FALSE & is.null(chooseYourGenes)) {
        warning(paste0("You don't have selected to plot gene names. 
                To plot only expreme points use 'extremePoints = TRUE' 
                and to plot all the gene names use 'allNames = TRUE'.
                You can find the plot in ", getwd()))
      }

      jpeg(outFile, width = 200, height = 200, units="mm", res=300)
      plot(pc, wmZ, xlim = c(0,1), xlab="P", ylab="z")

      ######################################## RECT = XLEFT, YBOTTOM, XRIGHT, YTOP
      ####################################
      # horizontal cut : max z-score = author cut : author max z-score
      ### NON-HUBS: z < 2.5 - proportional to  2.5*2.5/8 = 0.78125
      wmZ[is.na(wmZ)] <- 0
      maxZscore <- max(wmZ)
      hCut <- (maxZscore * 2.5) / 8
      minWmz <- min(wmZ-1)

      # for R1: ULTRA PERIPHERAL NODES (p<=0.05) = nodes with all their links within their module
      rect(-0.04, minWmz, 0.05, hCut, col="ivory4", border = NA)
      # for R2: PERIPHERAL NODES (0.05 < p <= 0.62) = nodes with most links within their module
      rect(0.05, minWmz, 0.62, hCut, col="lightsalmon", border = NA)
      # for R3: NON-HUB CONNECTOR NODES (0.62 < p <= 0.8) = nodes with many links to other modules
      rect(0.62, minWmz, 0.8, hCut, col = "darkseagreen", border = NA)
      # for R4: NON-HUB KINLESS NODES (p > 0.8) = nodes with links homogeneously distributed among all modules
      rect(0.8, minWmz, 1.3, hCut, col = "slateblue2", border = NA)
      ### HUBS : z > 2.5
      # for R5: PROVINCIAL HUBS (p <= 0.3) = hub nodes with the vast majority of links within their module
      rect(-0.3, hCut, 0.3, maxZscore+0.2, col ="lightyellow", border = NA)
      # for R6: CONNECTOR HUBS (0.3 < p <= 0.75) = hubs with many links to most of the other modules
      rect(0.3, hCut, 0.75, maxZscore+0.2, col = "mistyrose", border = NA)
      # for R7: KINLESS HUBS (p > 0.75) = hubs with links homogeneously distributed among all modules
      rect(0.75, hCut, 1.3, maxZscore+0.2, col = "gray90", border = NA)

      points(pc, wmZ, pch=16)

      text(c(0.03, 0.60, 0.78, 1.02, 0.28, 0.73, 1.02), 
           c(hCut, hCut, hCut, hCut, maxZscore+0.2, 
             maxZscore+0.2, maxZscore+0.2), paste0("R", 1:7), 
           cex=0.6, font=2, pos=1)

      if(any(vLabl!="")){
        plotrix::thigmophobe.labels(pc, wmZ, vLabl[match(names(pc), 
                                    V(topNetwork)$name)], cex=0.8)
      }

      dev.off()
    }
  }

  return(df)

}
