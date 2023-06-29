#' Symmetric Normalization of Adjancency Matrix
#' @details This function applies the following normalization: 
#' a_ij' = a_ij / sqrt(d_i d_j), where d_i is the degree of vertex i
#' @export
#' @param A adjacency matrix
#' @return normalized adjacency matrix
#' @usage normalize_adj_mat(A)
#' @examples 
#' \dontrun{normalize_adj_mat(A)}
#' 
normalize_adj_mat <- function(A=NULL){

  #dii = degree(i)
  #aij' = aij / sqrt(dii * dij)

  ArowSums <- sqrt(Matrix::rowSums(A))
  ArowSumsMat <- Matrix::Matrix(rep(ArowSums, dim(A)[1]), dim(A)[1], dim(A)[2])
  AcolSumsMat <- Matrix::Matrix(rep(ArowSums, dim(A)[1]), dim(A)[1], dim(A)[2], 
                        byrow=TRUE)
  Anorm <- A / ArowSumsMat / AcolSumsMat
  rownames(Anorm) <- rownames(A)
  colnames(Anorm) <- colnames(A)

  return(Anorm)
}
