#' Calculation of omega i (internal)
#' 
#' @description Calculation of omega i (internal)
#' @param idx index
#' @param Ai matrix
#' @param dSprod matrix with products of delta S
#' @return n

calc_omega_i <- function(idx=NULL, Ai=NULL, dSprod=NULL){

    AiIdx <- match(names(idx), rownames(Ai)) #order the matrix as idx
    AiPerm <- Ai[AiIdx, AiIdx]

    omegaVecti <- dSprod[1:length(idx), 1:length(idx)] * AiPerm
    omegaVecti <- sum(omegaVecti) / 2

    return(omegaVecti)

}
