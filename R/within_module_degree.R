#' Within-module degree
#'
#' @import igraph
#' @export

within_module_degree <-
  function (g = NULL,
            memb = NULL,
            A = NULL,
            useP = FALSE) {
    stopifnot(is_igraph(g))
    N <- max(memb)
    if (is.null(A))
      A <- as_adj(g, sparse = FALSE, names = FALSE)
    if (useP) {
      z <- rep.int(0, dim(A)[1L])
      z <- matrix(data = 0,
                  nrow = dim(A)[1L],
                  ncol = 3)
      for (S in seq_len(N)) {
        #print(S)
        deg <- A[memb == S,] %*% (memb == S)
        x <- unlist(lapply(deg, function(k) {
          # cat( sum(deg >= k), "\n")
          # cat(k, "\n", deg, "\n")
          p <- sum(deg >= k) / length(deg)
          return(p)
        }))
        z[memb == S, 2] <- x
        z[memb == S, 3] <- deg
      }
      z[, 1] <- -log10(z[, 2])
      colnames(z) <- c("wmd_score", "p", "wmd")
      #z[is.infinite(z)] <- 4#max(z[is.finite(z)]) + 0.1
      
    } else {
      z <- matrix(data = 0,
                  nrow = dim(A)[1L],
                  ncol = 4)
      Ki <- rep.int(0, dim(A)[1L])
      Ksi <- sigKsi <- rep.int(0, N)
      for (S in seq_len(N)) {
        x <- A[memb == S,] %*% (memb == S)
        Ki[memb == S] <- x
        Ksi[S] <- mean(x)
        sigKsi[S] <- sd(x)
      }
      z[, 1] <- (Ki - Ksi[memb]) / sigKsi[memb]
      z[is.infinite(z[, 1]), 1] <- 0
      z[, 2] <- Ki
      z[, 3] <- Ksi[memb]
      z[, 4] <- sigKsi[memb]
      colnames(z) <-
        c("wmd_score", "wmd", "mean_wmd", "var_wmd")
    }
    return(z)
    
  }
