\title{
Robust estimation of the concentration parameter kappa of the von Mises distribution from robust measures of dispersion
}\name{EstimateKappa}\alias{EstimateKappa}\keyword{ circular }\description{
Robust estimation of the concentration parameter of the von Mises distribution from robust measures of dispersion (Demni et al, 2026)
}\usage{
EstimateKappa(X, stat = c("CLMS","CMAD", "CLTS"), results_grid = NULL)
}\arguments{
  \item{X}{numeric; a vector of class \code{circular}}
  \item{stat}{character; one of the available estimators \code{c("CLMS","CMAD", "CLTS")}, respectively Circular Median Absolute Deviation, Circular Least Median Spread, Circular Least Trimmed Standard deviation.}
  \item{results_grid}{a dataframe of a custom grid of kappa used for interpolation. If NULL (default), the grid is handled automatically.}
}\value{
\item{CSD}{the estimated circular standard deviation based on the robust measures of dispersion}
\item{kappa}{the estimated kappa based on the estimated robust circular standard deviation}
}\author{Davide Buttarazzi, Houyem Demni
}\references{
Demni et al. (2026), Robust measures of dispersion for circular data with an anomaly detection rule. Under review.

%% ~put references to the literature/web site here ~
}\seealso{CircularViolinPlot}\examples{
#install.packages("circular")
require(circular)
data1 <- rvonmises(50, circular(pi/2), 5)
EstimateKappa(data1)
}
