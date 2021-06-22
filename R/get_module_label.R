#' Get module label
#'
#' @param xy layout coordinates
#' @param f factor describing modules
#' @param labs name vector of labels. Names  must be equal to unique values of f.
#' @export

get_module_label <- function(xy=NULL, f=NULL, labs=NULL){

  if(length(labs) != length(unique(f))){
    stop("ERROR: number of labels must be equal to the number of unique elements of f\n")
  }

  #create a unique df
  ans <- cbind(data.frame(id=1:nrow(xy), xy), f=f, stringsAsFactors = F)
  colnames(ans)[2:3] <- c("xx", "yy")

  ans <- merge(ans, labs, by.x="f", by.y=0, sort=F)
  colnames(ans)[ncol(ans)] <- "lab"
  ans$lab <- as.character(ans$lab)

  #find "centroids"
  xy_centers <- cbind(tapply(ans$xx, ans$f, mean), tapply(ans$yy, ans$f, mean))

  #find vertices close to centroids
  ans <- merge(ans, xy_centers, by.x="f", by.y=0)
  ans$dx <- (ans$V1 - ans$xx)^2
  ans$dy <- (ans$V2 - ans$yy)^2
  ans$d <- apply(ans[, c("dx", "dy")], 1, function(x) sum(x))
  ans_is.min <- tapply(ans$d, ans$f, min)

  ans <- merge(ans, ans_is.min, by.x="f", by.y=0)
  colnames(ans)[ncol(ans)] <- "isMin"

  ans$isMin <- ans$isMin  == ans$d
  ans$f_lab <- apply(ans[, c("f", "lab")], 1, paste0, collapse="_")

  #keep labels onyl for centroids
  ans$lab[!ans$isMin] <- ans$f_lab[!ans$isMin] <- ""
  ans$lab[duplicated(ans$f_lab)] <- ""

  return(ans[order(ans$id), c("f", "xx", "yy", "lab")])

}
