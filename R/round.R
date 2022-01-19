#' Scientific approximation
#'
#' @description scientific approximation
#' @param x numeric
#' @param y number of digits
#' @return truncated output
#' @examples 
#' \dontrun{round(x,y)}
#' @usage round(x,y)
#' @export

round <- function(x, y=0){
    return(trunc(x*10^y + 0.5)/(10^y))
}
