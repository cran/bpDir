\title{
Compute the robust measure of circular spread: the Circular Median Absolute Deviation (CMAD)
}\name{CircularMAD}\alias{CircularMAD}\keyword{ circular }\description{
Compute the Circular Median Absolute Deviation (Demni et al, 2026)
}\usage{
CircularMAD(X, h = FALSE)
}\arguments{
  \item{X}{numeric; a vector of class \code{circular} in \code{units=c("radians")}}
  \item{h}{logical; optionally produce histogram of the distances between the median and X, set to FALSE by default.}
}\value{
Circular Median Absolute Deviation (CMAD)
}\author{Davide Buttarazzi, Houyem Demni
}\references{
Demni et al. (2026), Robust measures of dispersion for circular data with an anomaly detection rule. Under review.

%% ~put references to the literature/web site here ~
}\seealso{CircularViolinPlot}\examples{
#install.packages("circular")
require(circular)
data1 <- rvonmises(50, circular(pi/2), 5)
CircularMAD(data1)
}
