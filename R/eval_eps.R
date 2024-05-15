#' Evaluate the effect of epsilon on NSI
#'
#' @description Evaluate the effect of epsilon on NSI
#' @param X0 initial condition matrix
#' @param Xs network diffusion steady state; required if W not provided
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, see normalize_adj_mat() function; required if Xs is not provided
#' @param eps optional vector with epsilon values
#' @param eps.n number of eps points for every order of magnitude
#' @param na number of orders of magnitudes above the one detected in X0 to set max eps
#' @param nb number of orders of magnitudes below the one detected in X0 to set min eps
#' @param G igraph object
#' @param top top genes to consider
#' @importFrom igraph graph_from_adjacency_matrix
#' @return list of data frames
#' @export
#'

eval_eps <- function(X0 = NULL, Xs = NULL, W=NULL, eps = NULL, G = NULL, top = NULL, eps.n=4, na=2, nb=2) {
    
    stopifnot(length(Xs)>0 | length(W)>0, length(G)>0 | length(W)>0, length(X0)>0)
    
    if (is.null(top)) {
        top <- min(300, nrow(X0))
    }
    
    if(is.null(G)){
        G <- graph_from_adjacency_matrix(sign(W), mode = "undirected")
    }
    
    if(is.null(eps)){
        cat("Generating epsilon based on X0 values...\n")
        l10maxX0 <- max(apply(X0, 2, function(x) max(log10(x)) + na))
        l10minX0 <- min(apply(X0, 2, function(x) min(log10(x[x>0])) - nb))
        oom <- round(l10maxX0 - l10minX0)
        cat("# orders of magnitude: ", oom, "\n")
        nval<-eps.n*round(l10maxX0 - l10minX0)
        cat("# values considered: ", nval, "\n")
        eps <- 10^seq(l10minX0, l10maxX0, length.out=nval)
    }
    cat("eps:", eps, "\n")
    
    
    if(is.null(Xs)){
        cat("Performing ND...\n")
        Xs <- ND(X0 = X0, W = W)
    }
    
    nsiEps <- vector("list", length(eps))
    csX0 <- nsiEps
    csXs <- nsiEps
    for (i in 1:length(eps)) {
        
        tempEps <- matrix( eps[i], ncol = ncol(X0), nrow = nrow(X0))
        
        #NSI
        nsiEps[[i]] <- Xs / (X0 + tempEps)
        
        # order by current NSI
        nsiOrder <- apply(nsiEps[[i]], 2, function(x) order(-x)[1:top])
        
        #omega on X0
        csX0[[i]] <- sapply(1:ncol(X0), function(j) omega(G, array(X0[, j], dimnames = list(rownames(X0)))[nsiOrder[, j]]))
        
        #omega on Xs
        csXs[[i]] <- sapply(1:ncol(Xs), function(j) omega(G, array(Xs[, j], dimnames = list(rownames(Xs)))[nsiOrder[, j]]))
        
    }
    
    aucX0 <- do.call(rbind, lapply(csX0, function(i_mat) colSums(i_mat)))
    aucXs <- do.call(rbind, lapply(csXs, function(i_mat) colSums(i_mat)))
    
    aucX0  <-  apply(aucX0, 2, function(x) x / max(x))
    aucXs  <-  apply(aucXs, 2, function(x) x / max(x))
    
    #optimal epsilon values
    opt_eps <- matrix(eps[sapply(1:ncol(X0), function(i) which.max(aucX0[, i] + aucXs[, i]))], ncol=ncol(X0))
                          
    cat("optimal eps:", opt_eps, "\n")
    
    out <-
        list(
            S = nsiEps,
            csX0 = csX0,
            csXs = csXs,
            aucX0 = aucX0,
            aucXs = aucXs,
            eps = eps,
            opt_eps = opt_eps
        )
    
    return(out)
    
}
