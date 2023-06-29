#' Assign communities
#' 
#' 
#' @export
#' 

assign_communities <- function(g=NULL, commList=NULL){
  
  V(g)$comm_id <- commList$membership[match(V(g)$name, commList$names)]
  V(g)$comm_id <- commList$membership
  print(table(V(g)$comm_id))
  
  return(g)
  
}