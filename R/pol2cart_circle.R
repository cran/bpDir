## function to convert from polar to cartesian coordinates
#' @keywords internal
#' @noRd

pol2cart_circle <- function(theta, r=1){
  for(i in 1:length(theta)){
    mat_cart <- matrix(NA, nrow=length(theta), ncol=2)
    #mat_cart <- c()
    x <- vector()
    y <- vector()
    x <- r*cos(theta)
    y <- r*sin(theta)
    mat_cart <- cbind(x, y)
  }
  return(mat_cart)
}
