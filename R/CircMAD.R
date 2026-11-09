#' @keywords internal
#' @noRd
CircMAD = function(X)
{
  #Compute the circular median for X
  SFMAngle=median.circular(as.circular(X))
  #Transform X and the circular median to cartesian coordinates
  X <- pol2cart_circle(X)
  SFM<-pol2cart_circle(SFMAngle)



  #Compute the arc distances of the Fisher median to each point in X
  arcdist <- acos(SFM%*%t(X))
  CMAD<-median(arcdist)
  return(CMAD)
}
