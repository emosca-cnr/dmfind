#' Linear map
#'
#' @description Linear map
#' @param x numeric vector
#' @param yMin min value of the new distribution
#' @param yMax max value of the new distribution
#' @export

linear_map <- function (x, yMin, yMax){
    m <- ((yMax - yMin)/(max(x) - min(x)))
    c <- yMax - m * max(x)
    y <- m * x + c
    return(y)
}
