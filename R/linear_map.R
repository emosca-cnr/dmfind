#' Linear map
#'
#' @description Linear map
#' @param x numeric vector
#' @param yMin min value of the new distribution
#' @param yMax max value of the new distribution
#' @return linear map y
#' @usage linear_map(x, yMin, yMax)
#' @examples linear_map(x, yMin, yMax)
#' @export

linear_map <- function (x, yMin, yMax){
    m <- ((yMax - yMin)/(max(x) - min(x)))
    c <- yMax - m * max(x)
    y <- m * x + c
    return(y)
}
