#' Enrichment Score
#'
#' @description function to calculate the enrichment score
#' @param idx vector of indices a subset of elements of x
#' @param x named vector, ranked list
#' @return enrichment score
#' @export
#'
es <- function(idx, x, le=F){

    #idx: indexes of a subset of elements of x
    #x: array of elements
    ## ES score
    #Phit(S,i) <- SUM_{i in idx, j=1..i}( r_j / Nr)
    #r: score
    #Nr: total score in x[idx]

    N <- length(x)
    Nh <- length(idx)
    hits <-  rep(0, length(x))
    hits[idx] <- x[idx]
    misses <- rep(1, length(x))
    misses[idx] <- 0
    hitsCumsum <- cumsum(abs(hits))
    Nr <- sum(abs(hits)) #equal to sum(abs(x[idx]))

    if(Nr==0){   #nothing to do
        es <- 0
        if(le){
            deviation <- 0
            tags <- 0
            tagsPerc <- 0
            listTop <- 0
            listTopPerc <- 0
            leadEdge <- 0
            leadEdgeSubset <- 0
        }
    }else{
        missesCumsum <- cumsum(misses)
        Phits <- hitsCumsum / Nr
        Pmiss <- missesCumsum / (N - Nh)
        deviation <- Phits - Pmiss
        wm <- which.max(abs(deviation))

        #es <- deviation[which.max(abs(deviation))]
        es <- deviation[wm]

        if(le){
            if(es >=0){
                tags <- sum(idx <= wm)
                listTop <- wm
                leadEdgeSubset <- paste0(names(x)[intersect(1:wm, idx)], 
                                           collapse = ";")
            }else{
                tags <- sum(idx > wm)
                listTop <- N - wm
                leadEdgeSubset <- paste0(names(x)[intersect(wm:N, idx)], 
                                           collapse = ";")
            }
            tagsPerc <- tags / length(idx)
            listTopPerc <- listTop / length(x)
        }
    }

    if(le){
        leadEdge <- tagsPerc * (1-listTopPerc) * N / (N - Nh)
        return(list(es=es, deviation=deviation, 
                    lea=data.frame(tags=tags, tagsPerc=tagsPerc, 
                    listTop=listTop, listTopPerc=listTopPerc, 
                    leadEdge=leadEdge, leadEdgeSubset=leadEdgeSubset, 
                    stringsAsFactors = F)))
    }else{
        return(es)
    }
}
