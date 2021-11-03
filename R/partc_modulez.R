#' pc_wmz Calculate participation coeffcient and within module degree z-score
#' @param top_network to do
#' @param all_names to do
#' @import brainGraph plotrix
#' @export

pc_wmz <- function(top_network, all_names = FALSE, extreme_points = FALSE, choose_your_genes = NULL, out_file="pc_wmz.jpg", do_plot=TRUE) {

  # partecipation coefficient
  pc <- brainGraph::part_coeff(top_network, V(top_network)$comm_id)

  # within module z-score
  wm_z <-brainGraph::within_module_deg_z_score(top_network, V(top_network)$comm_id)
  wm_z[is.na(wm_z)] <- 0


  # cbind pc and wm_z (x and y coordinates)
  df <- data.frame(name=names(pc), label=V(top_network)$label[match(names(pc), V(top_network)$name)], comm_id=V(top_network)$comm_id[match(names(pc), V(top_network)$name)], pc=pc, wm_z=wm_z, stringsAsFactors = F)

  max_zscore <- max(wm_z)
  h_cut <- (max_zscore * 2.5) / 8

  # look for coordinates
  # NON-HUBS : z < 2.5
  # for R1: ULTRA PERIPHERAL NODES (p<=0.05) = nodes with all their links within their module
  df$region[df$pc <= 0.05 & df$wm_z < h_cut] <- "R1"
  df$definition[df$region == "R1"] <- "ULTRA PERIPHERAL NODES"
  df$description[df$region == "R1"] <- "nodes with all their links within their module"

  # for R2: PERIPHERAL NODES (0.05 < p <= 0.62) = nodes with most links within their module
  df$region[df$pc > 0.05 & df$pc <= 0.62 & df$wm_z < h_cut] <- "R2"
  df$definition[df$region == "R2"] <- "PERIPHERAL NODES"
  df$description[df$region == "R2"] <- "nodes with most links within their module"

  # for R3: NON-HUB CONNECTOR NODES (0.62 < p <= 0.8) = nodes with many links to other modules
  df$region[df$pc > 0.62 & df$pc <= 0.8 & df$wm_z < h_cut] <- "R3"
  df$definition[df$region == "R3"] <- "NON-HUB CONNECTOR NODES"
  df$description[df$region == "R3"] <- "nodes with many links to other modules"

  # for R4: NON-HUB KINLESS NODES (p > 0.8) = nodes with links homogeneously distributed among all modules
  df$region[df$pc > 0.8 & df$wm_z < h_cut] <- "R4"
  df$definition[df$region == "R4"] <- "NON-HUB KINLESS NODES"
  df$description[df$region == "R4"] <- "nodes with links homogeneously distributed among all modules"

  ### HUBS : z > 2.5
  # for R5: PROVINCIAL HUBS (p <= 0.3) = hub nodes with the vast majority of links within their module
  df$region[df$pc <= 0.3 & df$wm_z > h_cut] <- "R5"
  df$definition[df$region == "R5"] <- "PROVINCIAL HUBS"
  df$description[df$region == "R5"] <- "hub nodes with the vast majority of links within their module"

  # for R6: CONNECTOR HUBS (0.3 < p <= 0.75) = hubs with many links to most of the other modules
  df$region[df$pc > 0.3 & df$pc <= 0.75 & df$wm_z > h_cut] <- "R6"
  df$definition[df$region == "R6"] <- "CONNECTOR HUBS"
  df$description[df$region == "R6"] <- "hubs with many links to most of the other modules"

  # for R7: KINLESS HUBS (p > 0.75) = hubs with links homogeneously distributed among all modules
  df$region[df$pc > 0.75 & df$wm_z > h_cut] <- "R7"
  df$definition[df$region == "R7"] <- "KINLESS HUBS"
  df$description[df$region == "R7"] <- "hubs with links homogeneously distributed among all modules"

  ########################################
  ########### check parameters ###########
  ########################################

  if(do_plot){
    if (all_names == TRUE & extreme_points == TRUE | all_names == TRUE & extreme_points == TRUE & !is.null(choose_your_genes)) {
      stop("Impossibile to satisfy all the conditions. Choose different parameters.")
    } else {

      ########################################
      ################# PLOT #################
      ########################################

      v_labl <- ""

      if (all_names == TRUE) {
        warning("You have selected to plot all gene names. To plot only expreme points use 'extreme_points = TRUE")
        v_labl <- V(top_network)$label
      }
      if (extreme_points == TRUE) {
        warning("You have selected to plot extreme gene names. To plot all the gene names use 'all_names = TRUE'")

        name_max_zscore <- V(top_network)$label[match(names(pc)[which.max(wm_z)], V(top_network)$name)]

        # gene with max pc
        name_max_pc <- V(top_network)$label[match(names(pc)[which.max(pc)], V(top_network)$name)]

        # genes with pc = 0
        V(top_network)$label[match(names(pc)[pc == 0 & which.max(wm_z)], V(top_network)$name)]

        # gene with min z score
        name_min_zscore <- V(top_network)$label[match(names(pc)[which.min(wm_z)], V(top_network)$name)]

        # gene with min pc
        name_min_pc <- V(top_network)$label[match(names(pc)[which.min(pc)], V(top_network)$name)]

        v_labl <- V(top_network)$label
        v_labl[!v_labl %in% c(name_max_zscore, name_max_pc, name_min_zscore, name_min_pc)] <- ""
      }

      if (!is.null(choose_your_genes)) {

        warning("You have selected to plot your gene names. To plot all the gene names use 'all_names = TRUE'")
        v_labl <- V(top_network)$label
        v_labl[!v_labl %in% choose_your_genes] <- ""

      }

      #V(top_network)$v_labl <- v_labl

      # PLOT ALL GENE NAMES
      if (all_names == FALSE & extreme_points == FALSE & is.null(choose_your_genes)) {
        warning(paste0("You don't have selected to plot gene names. To plot only expreme points use
                     'extreme_points = TRUE' and to plot all the gene names use 'all_names = TRUE'.
                     You can find the plot in ", getwd()))
      }


      jpeg(out_file, width = 200, height = 200, units="mm", res=300)

      plot(pc, wm_z, xlim = c(0,1), xlab="P", ylab="z")

      ######################################## RECT = XLEFT, YBOTTOM, XRIGHT, YTOP
      ####################################
      # horizontal cut : max z-score = author cut : author max z-score
      ### NON-HUBS: z < 2.5 - proportional to  2.5*2.5/8 = 0.78125
      wm_z[is.na(wm_z)] <- 0
      max_zscore <- max(wm_z)
      h_cut <- (max_zscore * 2.5) / 8

      min_wmz <- min(wm_z-1)

      # for R1: ULTRA PERIPHERAL NODES (p<=0.05) = nodes with all their links within their module
      rect(-0.04, min_wmz, 0.05, h_cut, col="ivory4", border = NA)
      # for R2: PERIPHERAL NODES (0.05 < p <= 0.62) = nodes with most links within their module
      rect(0.05, min_wmz, 0.62, h_cut, col="lightsalmon", border = NA)
      # for R3: NON-HUB CONNECTOR NODES (0.62 < p <= 0.8) = nodes with many links to other modules
      rect(0.62, min_wmz, 0.8, h_cut, col = "darkseagreen", border = NA)
      # for R4: NON-HUB KINLESS NODES (p > 0.8) = nodes with links homogeneously distributed among all modules
      rect(0.8, min_wmz, 1.3, h_cut, col = "slateblue2", border = NA)
      ### HUBS : z > 2.5
      # for R5: PROVINCIAL HUBS (p <= 0.3) = hub nodes with the vast majority of links within their module
      rect(-0.3, h_cut, 0.3, max_zscore+0.2, col ="lightyellow", border = NA)
      # for R6: CONNECTOR HUBS (0.3 < p <= 0.75) = hubs with many links to most of the other modules
      rect(0.3, h_cut, 0.75, max_zscore+0.2, col = "mistyrose", border = NA)
      # for R7: KINLESS HUBS (p > 0.75) = hubs with links homogeneously distributed among all modules
      rect(0.75, h_cut, 1.3, max_zscore+0.2, col = "gray90", border = NA)

      points(pc, wm_z, pch=16)

      text(c(0.03, 0.60, 0.78, 1.02, 0.28, 0.73, 1.02), c(h_cut, h_cut, h_cut, h_cut, max_zscore+0.2, max_zscore+0.2, max_zscore+0.2), paste0("R", 1:7), cex=0.6, font=2, pos=1)

      if(any(v_labl!="")){
        plotrix::thigmophobe.labels(pc, wm_z, v_labl[match(names(pc), V(top_network)$name)], cex=0.8)
      }

      dev.off()


    }
  }

  return(df)
  
  # changes

}
