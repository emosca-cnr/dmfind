#' Calculation of omega i (internal)
#' 
#' @description Calculation of omega i (internal)
#' @param idx index
#' @param Ai matrix
#' @param dSprod matrix with products of delta S
#' @param norm normalize by number of links (l), number or vertices (v) or do not normalize (n)

calc_omega_i <- function(idx, Ai, dSprod, norm=c("n", "l", "v")){

    norm <- match.arg(norm, c("n", "l", "v"))

    Ai.idx <- match(names(idx), rownames(Ai)) #order the matrix as idx
    Ai_perm <- Ai[Ai.idx, Ai.idx]

    omega_vect_i <- dSprod[1:length(idx), 1:length(idx)] * Ai_perm
    omega_vect_i <- sum(omega_vect_i) / 2

    if(norm=="l" & (omega_vect_i > 0)){
        Ai_perm[upper.tri(Ai_perm)] <- 0
        nE <- sum(Ai_perm)
        omega_vect_i <- omega_vect_i / nE
    }
    if(norm=="v"){
        omega_vect <- omega_vect / length(omega_vect)
    }

    return(omega_vect_i)

}
