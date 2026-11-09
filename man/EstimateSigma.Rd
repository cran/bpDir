\title{
Robust estimation of the scale parameter sigma of the Wrapped-Normal distribution from the robust measures of dispersion
}\name{EstimateSigma}\alias{EstimateSigma}\keyword{ circular }\description{
Robust estimation of the scale parameter sigma of the Wrapped-Normal distribution from the robust measures of dispersion (Demni et al, 2026)
}\usage{
EstimateSigma(X, stat = c("CLMS", "CMAD", "CLTS"), results_grid = NULL)
}\arguments{
  \item{X}{numeric; a vector of class \code{circular}}
  \item{stat}{character; one of the available estimators \code{c("CLMS", "CMAD", "CLTS")}, respectively Circular Median Absolute deviation, Circular Least Median Spread, Circular Least Trimmed Standard deviation.}
  \item{results_grid}{a dataframe of a custom grid of kappa used for interpolation. If NULL (default), the grid is handled automatically.}
}\value{
\item{CSD}{the estimated circular standard deviation based on the robust measures of dispersion}
\item{sigma}{the estimated sigma based on the estimated robust circular standard deviation}
}\author{Davide Buttarazzi, Houyem Demni
}\references{
Demni et al. (2026), Robust measures of dispersion for circular data with an anomaly detection rule. Under review.

%% ~put references to the literature/web site here ~
}\seealso{CircularViolinPlot}\examples{
#install.packages("circular")
require(circular)
data1 <- rvonmises(25, circular(pi/2), 5)
g=generate_results_grid_sigma(sigma_seq = seq(0.1, 1.5, by = 0.1))
EstimateSigma(data1, results_grid = g)
}
