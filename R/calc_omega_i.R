#' Calculation of omega i (internal)
#' 
#' @description Calculation of omega i (internal)
#' @param idx index
#' @param Ai matrix
#' @param dSprod matrix with products of delta S
#' @param norm normalize by number of links (l), number or 
#' vertices (v) or do not normalize (n)
#' @return n

calc_omega_i <- function(idx, Ai, dSprod, norm=c("n", "l", "v")){

    norm <- match.arg(norm, c("n", "l", "v"))

    AiIdx <- match(names(idx), rownames(Ai)) #order the matrix as idx
    AiPerm <- Ai[AiIdx, AiIdx]

    omegaVecti <- dSprod[1:length(idx), 1:length(idx)] * AiPerm
    omegaVecti <- sum(omegaVecti) / 2

    if(norm=="l" & (omegaVecti > 0)){
        AiPerm[upper.tri(AiPerm)] <- 0
        nE <- sum(AiPerm)
        omegaVecti <- omegaVecti / nE
    }
    if(norm=="v"){
        omegaVect <- omegaVect / length(omegaVect)
    }

    return(omegaVecti)

}
