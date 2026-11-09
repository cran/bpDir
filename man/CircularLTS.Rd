\title{
Compute the robust measure of circular spread: the Circular Least Trimmed Standard deviation (CLTS)
}\name{CircularLTS}\alias{CircularLTS}\keyword{ circular }\description{
Compute the Circular Least Trimmed Standard deviation (Demni et al, 2026)
}\usage{
CircularLTS(X)
}\arguments{
  \item{X}{numeric; a vector of class \code{circular} in \code{units=c("radians")}}
}\value{
Circular Least Trimmed Standard deviation (CLTS)
}\author{Davide Buttarazzi, Houyem Demni
}\references{
Demni et al. (2026), Robust measures of dispersion for circular data with an anomaly detection rule. Under review.

%% ~put references to the literature/web site here ~
}\seealso{CircularViolinPlot}\examples{
#install.packages("circular")
require(circular)
data1 <- rvonmises(50, circular(pi/2), 5)
CircularLTS(data1)
}
