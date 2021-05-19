#' Local Monte Carlo p values
#' @param G gene x gene undirected interaction graph
#' @param ranked_vector scores vector; it must have the same names and size of the vertices of G
#' @param k number of permutations
#' @param mc.cores number of cores
#' @param null_model null model type
#' @export
#' @import parallel


NR <- function(G, ranked_vector, k=99, mc.cores=1, null_model="A", norm="n") {

  null_model <- match.arg(null_model)

  #real values
  ranked_vector_norm <- ranked_vector / max(ranked_vector)

  omega_vect <- omega(G, ranked_vector_norm, norm = norm)

  #permutations of ranked_vector_norm: named vector of indices
  idx.perm <- lapply(1:k, function(x) array(sample(1:length(ranked_vector_norm)), dimnames = list(names(ranked_vector_norm))))

  if(mc.cores > 1){
    res <- parallel::mclapply(idx.perm, omega_perm, G=G, dS=ranked_vector_norm, mc.cores = mc.cores, null_model=null_model, norm=norm)
  }else{
    res <- lapply(idx.perm, omega_perm, G=G, dS=ranked_vector_norm, null_model=null_model, norm=norm)
  }

  #n-by-k matrix
  res <- t(do.call(rbind, res))

  #LMC p values
  out <- apply(res, 2, function(x) sign(x >= omega_vect))
  out <- rowSums(out)
  out <- (out + 1) / (k+1)

  return(list(NR_summary=data.frame(id=names(out), rank=1:length(ranked_vector_norm), ranking_score=ranked_vector_norm, omega=omega_vect, p=out, stringsAsFactors = FALSE), omega_perm=res))

}


