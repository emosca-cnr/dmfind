#' Calculate ccf
#' @export

calculate_ccf <- function(G=NULL, ranked_vector_names=NULL){

    n <- length(ranked_vector_names)
    ranked_vector_names <- 
        ranked_vector_names[1:min(length(ranked_vector_names), 1000)]
    n <- length(ranked_vector_names)

    ans <- sapply(1:n, function(i) 
        CCF(igraph::induced_subgraph(g, 
        V(g)$name[V(g)$name %in% ranked_vector_names[1:i]])))

    return(ans)

}
