#' calc_Mcp
#' 
#' @description Help function
#' @param perms list of permuted statistics
#' @param real vector of real statistics
#' @return n
#'
#' 
.calc_MCp <- function(perms, real){

    check_extremes <- function(xxx, yyy){
        #xxx real
        #yyy permuted
        if(!is.numeric(yyy)){
          print(yyy)
          stop('not numeric\n')
        }
    
        out <- list(
            g = which(yyy >= xxx),
            l = which(yyy <= xxx),
            abs= which(abs(yyy) >= abs(xxx)))
    
        return(out)
    }

    temp <- lapply(perms, function(xx) check_extremes(real, xx))
    pPerms <- list(
        g = names(real)[unlist(lapply(temp, function(xx) xx$g))],
        l = names(real)[unlist(lapply(temp, function(xx) xx$l))],
        abs = names(real)[unlist(lapply(temp, function(xx) xx$abs))])

    nperms <- length(perms) + 1

    zeroCases <- vector('list', 3)
    names(zeroCases) <- names(pPerms)
    for(i in 1:3){
        zeroCases[[i]] <- names(real)[!(names(real) %in% pPerms[[i]])]
        temp <- rep(1, length(zeroCases[[i]]))
        names(temp) <- zeroCases[[i]]
        zeroCases[[i]] <- temp / (nperms)

        pPerms[[i]] <- (table(pPerms[[i]]) + 1) / (nperms)
        pPerms[[i]] <- c(pPerms[[i]], zeroCases[[i]])
    }
    rm(temp)

    out <- merge(pPerms$g, pPerms$l, by=0, sort=FALSE)
    colnames(out)[colnames(out)=='x'] <- 'p_g'
    colnames(out)[colnames(out)=='y'] <- 'p_l'

    out <- merge(out, pPerms$abs, by.x=1, by.y=0, sort=FALSE)
    colnames(out)[colnames(out)=='y'] <- 'p_abs'

    rownames(out) <- out$Row.names
    out$Row.names <- NULL

    return(out[match(names(real), rownames(out)),])

}
