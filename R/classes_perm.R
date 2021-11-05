#' Generate N permutations of a binary vector
#' 
#' @description Generate N permutations of a binary vector
#' @param z binary vector with values in \{1, 2\}
#' @param N number of permutations
#' @return list of N permutations
#' @export
#' 
classes_perm <- function(z, N){

    N2 <- length(which(z == 2))

    z.perms <- lapply(1:N, function(x) rep(1, length(z)))
    z.perms <- lapply(z.perms, function(x) 
        replace(x, sample(1:length(x), N2), 2))

    return(z.perms)

}
