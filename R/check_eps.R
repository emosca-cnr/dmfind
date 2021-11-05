#' Search for the optimal epsilon
#' 
#' @description Search for optimal epsilon
#' @param X0 1-column matrix
#' @param Xs 1-column matrix (smoothed)
#' @param x optional vector with top names to use, in place of X0, in the comparison with S
#' @param eps values to check
#' @param top number of top genes on which the comparison has to be done
#' @param sort.X0 whether to sort (TRUE) or not (FALSE) X0 in decreasing order
#' @export
#'
check_eps <- function (X0, Xs, x=NA, eps = c(0.1, 0.25, 0.5, 1), 
                       top = 250, sort.X0=FALSE){
    #res <- data.frame(eps = eps, oi = 0, stringsAsFactors = F)
    S_eps <- lapply(eps, function(xx) nsi(X0, Xs, eps = xx))
    S_eps <- lapply(S_eps, function(xx) xx[order(-xx[, 1])[1:top], , drop = F])

    #sort X0
    X0_sorted <- X0
    if(sort.X0){
        X0_sorted <- X0_sorted[order(-X0_sorted[, 1]), , drop = F]
    }
    X0_sorted <- X0_sorted[1:top, , drop=F]
    N_pos_X0 <- sum(X0_sorted>0)
    #intersection
    res <- vector("list", length(eps))
    for(i in 1:length(eps)){
        res[[i]] <- unlist(lapply(1:top, function(xx) 
            sum(rownames(X0_sorted)[X0_sorted>0] %in% 
                    rownames(S_eps[[i]])[1:xx])/min(xx, N_pos_X0)))
    }
    res <-do.call(cbind, res)
    colnames(res) <- paste0("eps", eps)
    return(res)
}
