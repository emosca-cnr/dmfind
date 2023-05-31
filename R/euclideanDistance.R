euclideanDistance <- function(Xa, Xb, Ya, Yb) {
    
    euDist <- sqrt((Xb-Xa)^2+(Yb-Ya)^2)
    return(euDist)
    
}