#' Gene Set Enrichment Analysis
#' 
#' @description Gene Set Enrichment Analysis
#' @param rl numeric matrix of genes-by-ranking criteria; each column contains 
#' numeric values; rownames are mandatory
#' @param gsl named list of gene sets
#' @param k integer, number of permutations
#' @param ordMode ordering mode: -1 -> descending; 1 ascending; must be 
#' of length equal to ncol(rl)
#' @param mcCoresPath number of cores to use for parallel calculation of gene 
#' set lists; the total number of cpu used will be mcCoresPath x mcCoresPerm
#' @param mcCoresPerm number of cores to use for parallel calculation of ranked 
#' list permutations; the total number of cpu used will be mcCoresPath x mcCoresPerm
#' @import BiocParallel
#' @return data.frame with es, nes, p-value, adjusted p-value and FDR q-value
#' @examples 
#' \dontrun{gsea(rl, gsl, k=100, ordMode=-1, mcCoresPath=1, mcCoresPerm=1)}
#' @usage gsea(rl, gsl, k=100, ordMode=-1, mcCoresPath=1, mcCoresPerm=1)
#' @export

gsea <- function(rl, gsl, k=100, ordMode=-1, 
                 mcCoresPath=1, mcCoresPerm=1){

    #cheks
    if(!is.matrix(rl) | !is.numeric(rl)){
        stop("rl must be a numeric matrix")
    }

    if(length(ordMode) != ncol(rl)){
        stop("length(ordMode) must be equal to ncol(rl)")
    }

    #create the list of ranked vectors
    rll <- vector('list', ncol(rl))
    names(rll) <- colnames(rl)
    for(i in 1:length(rll)){
        rll[[i]] <- sort(array(rl[, i], 
                         dimnames = list(rownames(rl))), 
                         decreasing = ordMode[i]==-1)
    }

    #permutation of gene ids
    cat('generating', k, 'permutations\n')
    xPerm <- lapply(1:k, function(x) sample(rownames(rl), nrow(rl)))

    #real es
    print("ES...")
    #realEs <- lapply(gsl, function(x) lapply(rll, function(y) 
    #es(which(names(y) %in% x), y)))
    #realEs <- do.call(rbind, realEs)

    realEsData <- lapply(gsl, function(x) 
        lapply(rll, function(y) es(which(names(y) %in% x), y, le=TRUE)))
    realEs <- do.call(rbind, 
                       lapply(realEsData, 
                              function(x) unlist(lapply(x, function(y) y$es))))


    leadingEdge <- vector("list", length(rll))
    names(leadingEdge) <- names(rll)
    for(i in 1:length(leadingEdge)){
        leadingEdge[[i]] <- 
            do.call(rbind, lapply(realEsData, function(x) x[[i]]$lea))
    }

    #permutations
    print("calculating permutations")
    if(mcCoresPath==1){
        if(mcCoresPerm == 1){
        #res <- lapply(gsl, function(x) unlist(lapply(xPerm, 
            #function(y) calc_gs_perm(rl, y, x))))
            res <- lapply(gsl, function(x) 
                do.call(rbind, lapply(xPerm, function(y) 
                    calc_gs_perm(rll, y, x))))
        }else{
            cat(k, "permutations on", mcCoresPerm, "cores\n")
            #res <- lapply(gsl, function(x) unlist(mclapply(xPerm, 
            #function(y) calc_gs_perm(rl, y, x), mc.cores=mcCoresPerm)))
            res <- lapply(gsl, function(x) 
                    do.call(rbind, BiocParallel::bplapply(xPerm, function(y) 
                        calc_gs_perm(rll, y, x), mc.cores=mcCoresPerm)))
        }
    }else{
        cat(length(gsl), "gene sets on", mcCoresPath, "cores\n")
        if(mcCoresPerm == 1){
            #res <- parallel::mclapply(gsl, function(x) 
            #unlist(lapply(xPerm, function(y) calc_gs_perm(rl, y, x))), 
            #mc.cores = mcCoresPath)
            res <- parallel::mclapply(gsl, function(x) 
                    do.call(rbind, lapply(xPerm, function(y) 
                    calc_gs_perm(rll, y, x))), mc.cores = mcCoresPath)
        }else{
            cat(k, "permutations on", mcCoresPerm, "cores\n")
            #res <- parallel::mclapply(gsl, function(x) 
            #unlist(parallel::mclapply(xPerm, function(y) 
            #calc_gs_perm(rl, y, x), mc.cores=mcCoresPerm)), 
            #mc.cores = mcCoresPath)
            res <- parallel::mclapply(gsl, function(x) 
                    do.call(rbind, BiocParallel::bplapply(xPerm, function(y) 
                    calc_gs_perm(rll, y, x), mc.cores=mcCoresPerm)), 
                    mc.cores = mcCoresPath)
        }
    }

    #list of pathwa-by-k matrices of permuted es
    #res <- do.call(rbind, res)
    #   if(nrow(res) != length(gsl) | ncol(res) != k){
    #     stop("not all the permutations returned a correct value\n")
    #   }
    #res <- cbind(realEs, res)

    temp <- vector('list', length(rll))
    names(temp) <- colnames(rl)
    for(i in 1:length(rll)){
        temp[[i]] <- cbind(realEs[, i], 
                           do.call(rbind, lapply(res, function(x) x[, i])))
    }
    res <- temp
    rm(temp)

    #statistics
    print("calculating statistics...")

    #the first column is the real value, and it is 
    #included in the calculation of p
    #if ES* > 0 -> p = # (ESp >= ES*) / (k+1)
    #if ES* < 0 -> p = # (ESp <= ES*) / (k+1)
    out <- res

    for(i in 1:length(rll)){
        pVal <- apply(res[[i]], 1, function(x) 
                    ifelse(x[1] >= 0, sum(x >= x[1]) / length(x), 
                    sum(x <= x[1]) / length(x)))
        #normalized ES
        #positive, negative
        means <- t(apply(res[[i]], 1, function(x) 
            c(mean(x[x>0]), abs(mean(x[x<0]))))) 
        nes <- res[[i]] / means[, 1]
        nesNeg <- res[[i]] / means[, 2]
        nes[res[[i]] < 0] <- nesNeg[res[[i]] < 0]
        rm(means, nesNeg)
        #calculate FDR
        allNes <- as.numeric(nes)
        nNesPos <- sum(nes>=0)
        nNesNeg <- sum(nes<=0)
        nRealNesPos <- sum(nes[,1] >= 0)
        nRealNesNeg <- sum(nes[,1] <= 0)
        if((nNesPos + nNesNeg) != (k*length(gsl) + nrow(nes))){
            warning('some nes value is equal to zero')
        }

        #FDR: NES* > 0: fdrq = #(all positive NESp >= NES*) / #(all positive NESp) / [ #(all NES* >= NES*) / (all positive NES*) ]
        #FDR: NES* < 0: fdrq = #(all negative NESp <= NES*) / #(all negative NESp) / [ #(all NES* <= NES*) / (all negative NES*) ]
        fdrq <- sapply(nes[, 1], function(x) 
                ifelse(x>0, sum(allNes >= x) / nNesPos, 
                sum(allNes <= x) / nNesNeg))
        fdrq <- fdrq / sapply(nes[, 1], function(x) 
                ifelse(x>0, sum(nes[, 1] >= x) / nRealNesPos, 
                sum(nes[, 1] <= x) / nRealNesNeg))
        rm(allNes)

        #out table
        out[[i]] <- data.frame(id=rownames(res[[i]]), es=res[[i]][, 1], 
                               pVal=pVal, 
                               adjpVal=stats::p.adjust(pVal, method='fdr'), 
                               nes=nes[, 1], FDRq=fdrq, stringsAsFactors=FALSE)

    }

    print("done")
    return(list(gsTable=out, leadingEdge=leadingEdge))

}
