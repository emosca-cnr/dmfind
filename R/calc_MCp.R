calc_MCp <- function(perms, real){

  #perms: list of permuted statistics
  #real: vector of real statistics

  check.extremes <- function(xxx, yyy){

    #xxx real
    #yyy permuted

    if(!is.numeric(yyy)){
      print(yyy)
      stop('not numeric\n')
    }

    out <- list(
      g = which(yyy >= xxx),
      l = which(yyy <= xxx),
      abs= which(abs(yyy) >= abs(xxx))
    )

    return(out)
  }

  temp <- lapply(perms, function(xx) check.extremes(real, xx))
  p.perms <- list(
    g = names(real)[unlist(lapply(temp, function(xx) xx$g))],
    l = names(real)[unlist(lapply(temp, function(xx) xx$l))],
    abs = names(real)[unlist(lapply(temp, function(xx) xx$abs))]
  )

  nperms <- length(perms) + 1

  zero.cases <- vector('list', 3)
  names(zero.cases) <- names(p.perms)
  for(i in 1:3){
    zero.cases[[i]] <- names(real)[!(names(real) %in% p.perms[[i]])]
    temp <- rep(1, length(zero.cases[[i]]))
    names(temp) <- zero.cases[[i]]
    zero.cases[[i]] <- temp / (nperms)

    p.perms[[i]] <- (table(p.perms[[i]]) + 1) / (nperms)
    p.perms[[i]] <- c(p.perms[[i]], zero.cases[[i]])
  }
  rm(temp)

  out <- merge(p.perms$g, p.perms$l, by=0, sort=FALSE)
  colnames(out)[colnames(out)=='x'] <- 'p_g'
  colnames(out)[colnames(out)=='y'] <- 'p_l'

  out <- merge(out, p.perms$abs, by.x=1, by.y=0, sort=FALSE)
  colnames(out)[colnames(out)=='y'] <- 'p_abs'

  rownames(out) <- out$Row.names
  out$Row.names <- NULL

  return(out[match(names(real), rownames(out)),])

}
