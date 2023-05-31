#' Generate N permutations of a binary vector
#' 
#' @description Generate N permutations of a binary vector
#' @param z binary vector with values in \{1, 2\}
#' @param N number of permutations
#' @return list of N permutations
#' @examples 
#' \dontrun{classes_perm(z,N)}
#' @usage classes_perm(z,N)
#' @export
#' 
classes_perm <- function(z, N){

    N2 <- length(which(z == 2))

    zPerms <- lapply(1:N, function(x) rep(1, length(z)))
    zPerms <- lapply(zPerms, function(x) 
        replace(x, sample(1:length(x), N2), 2))

    return(zPerms)

}
