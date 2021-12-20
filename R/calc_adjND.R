#' Calculation of permutation-adjusted network smoothing index
#' 
#' @param X0 input matrix
#' @param W symmetrically normalized adjacency matrix W = D^-1 A D^-1, see normalize_adj_mat function
#' @param k number of permutations
#' @param eps numeric value
#' @param mc.cores number of cores
#' @param mode whether to use network smoothing index (S) or network smoothing values (Xs)
#' @param ... additional parameteres of ND
#' @description Calculation of permutation-adjusted network smoothing index
#' @usage calc_adjND(X0, W, eps, k=99, mode=c("S", "Xs"), mc.cores=1, return.perm=FALSE)
#' @examples calc_adjND(X0, W, eps, k=99, mode=c("S", "Xs"), mc.cores=1, return.perm=FALSE)
#' @return \code{data.frame} with X0, Xs, S, p and Sp
#' @import BiocParallel
#' @export
#'
calc_adjND <- function(X0, W, eps=rep(1, ncol(X0)), k=99, mode=c("S", "Xs"), 
                       mc.cores=1, return.perm=FALSE, ...){

    mode <- match.arg(mode)

    allX0 <- c(list(X0), lapply(1:k, function(x) matrix(as.numeric(X0), 
                                ncol=ncol(X0), dimnames = 
                                list(sample(rownames(X0), nrow(X0))))))
    allX0 <- BiocParallel::bplapply(allX0, function(x) 
        x[match(rownames(W), rownames(x)), , drop=FALSE ])

    cat("network propagation\n")
    if(mc.cores==1){
        allXS <- lapply(allX0, function(x) ND(x, W, ...)$Xs)
    }else{
        allXS <- BiocParallel::bplapply(allX0, function(x) 
          ND(x, W, ...)$Xs, mc.cores=mc.cores)
    }


    if(mode=="S"){
        cat("calculation of Sp\n")
        #S
        all_S <- lapply(1:(k+1), function(i) nsi(allX0[[i]], 
                                                 allXS[[i]], eps=eps))
        Xs <- allXS[[1]]
    }else{
        all_S <- allXS
    }

    ## trash all permutations
    if(!return.perm){
        rm(allX0, allXS)
    }

    #p
    estP <- calc_p(all_S)
    S <- all_S[[1]]

    if(mode=="S"){
        if(return.perm){
            out <- list(Xs=Xs, eps=eps, S=S, p=estP, Sp= S * -log10(estP),
                        allX0=allX0, allXS=allXS)
        }else{
            out <- list(Xs=Xs, eps=eps, S=S, p=estP, Sp= S * -log10(estP))
        }
    }else{
        if(return.perm){
            out <- list(Xs=S, eps=eps, p=estP, Xsp= S * -log10(estP), 
                        allX0=allX0, allXS=allXS)
        }else{
            out <- list(Xs=S, eps=eps, p=estP, Xsp= S * -log10(estP))
        }
    }
    rm(all_S) ## trash all permutations

    return(out)
}
