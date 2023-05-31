#' cmp_pc_wmz Comparison between functional cartographies
#'
#' @description Comparison between functional cartographies
#' @param PcWmz_A dataframe 1 resulting from pc_wmz
#' @param PcWmz_B dataframe 2 resulting from pc_wmz
#' @return dataframe with deltaP, deltaZ, changeRegion, euclideanDistance, 
#' pc and wm_z from first and second network
#' @examples
#' \dontrun{cmp_pc_wmz(PcWmz_A, PcWmz_B)}
#' @usage cmp_pc_wmz(PcWmz_A, PcWmz_B)
#' @export

cmp_pc_wmz <- function(PcWmz_A, PcWmz_B) {
    
    cmpDf <- merge(firstPcWmz, secondPcWmz, by=c("name", "label"), all=T)
    cmpDf$euclideanDistance <- euclideanDistance(cmpDf$pc.x, cmpDf$pc.y, 
                                         cmpDf$wm_z.x, cmpDf$wm_z.y)
    cmpDf$changeRegion <- ifelse(cmpDf$region.x == cmpDf$region.y, "NO", "YES")
    
    colnames(cmpDf)[which(names(cmpDf) == "pc.x")] <- "pc_A"
    colnames(cmpDf)[which(names(cmpDf) == "wm_z.x")] <- "wm_z_A"
    colnames(cmpDf)[which(names(cmpDf) == "pc.y")] <- "pc_B"
    colnames(cmpDf)[which(names(cmpDf) == "wm_z.y")] <- "wm_z_B"
    colnames(cmpDf)[which(names(cmpDf) == "region.x")] <- "region_A"
    colnames(cmpDf)[which(names(cmpDf) == "region.y")] <- "region_B"
    colnames(cmpDf)[which(names(cmpDf) == "definition.x")] <- "definition_A"
    colnames(cmpDf)[which(names(cmpDf) == "definition.y")] <- "definition_B"
    colnames(cmpDf)[which(names(cmpDf) == "description.x")] <- "description_A"
    colnames(cmpDf)[which(names(cmpDf) == "description.y")] <- "description_B"
    colnames(cmpDf)[which(names(cmpDf) == "comm_id.x")] <- "comm_id_A"
    colnames(cmpDf)[which(names(cmpDf) == "comm_id.y")] <- "comm_id_A"
    
    cmpDf$deltaP <- cmpDf$pc_A-cmpDf$pc_B
    cmpDf$deltaZ <- cmpDf$wm_z_A-cmpDf$wm_z_B
    
    return(cmpDf)
    
}