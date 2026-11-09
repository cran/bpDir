######################### Ko estimator
#' @export
EstimateKappaKo <-function(X)
{
  CMAD<-CircMAD(X)
  kappa<- (qnorm(0.75)/(CMAD))^2
  return(kappa)
}
