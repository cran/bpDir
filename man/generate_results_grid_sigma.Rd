\name{generate_results_grid_sigma}
\alias{generate_results_grid_sigma}
\title{Generate a Calibration Grid for Estimating the Wrapped Normal Scale Parameter}
\description{
Generates a calibration grid relating robust measures of circular
dispersion to the circular standard deviation (CSD) and the scale
parameter (\eqn{\sigma}) of a wrapped normal distribution.

The grid is used by the estimation procedure to convert an observed
robust circular dispersion statistic into an estimate of CSD through
monotone spline interpolation, and subsequently into an estimate of
the wrapped normal scale parameter.
}
\usage{
generate_results_grid_sigma(
  sigma_seq = seq(0.1, 1.5, by = 0.0005)
)
}
\arguments{
  \item{sigma_seq}{
    A numeric vector of positive values of the wrapped normal scale
    parameter \eqn{\sigma} for which the calibration grid is generated.
    The default sequence ranges from 0.1 to 1.5 in increments of 0.0005.
  }
}
\value{
A data frame with the following columns:
\describe{
  \item{sigma}{Values of the wrapped normal scale parameter.}
  \item{rho}{The mean resultant length of the wrapped normal distribution,
  calculated as \eqn{\exp(-\sigma^2/2)}.}
  \item{CSD}{The circular standard deviation, calculated as
  \eqn{\sqrt{2(1-\rho)}}.}
  \item{CMAD}{The calibration values obtained using the 0.75 quantile
  of the centred wrapped normal distribution, for the CMAD-based
  estimation procedure.}
  \item{CLMS}{The calibration values obtained using the 0.75 quantile
  of the centred wrapped normal distribution, for the CLMS-based
  estimation procedure.}
  \item{CLTS}{The calibration values obtained using the CLTS-based estimation procedure.}
}
}
\details{
The calibration grid is constructed by evaluating the dispersion
statistics over a sequence of wrapped normal distributions with
different values of \eqn{\sigma}.

For a given observed dataset, the estimation procedure first computes
one of the following robust measures of circular dispersion:
\code{CircularMAD}, \code{CircularLMS}, or \code{CircularLTS},
depending on the selected statistic.

A precomputed calibration grid can be supplied to the estimation
function through its \code{results_grid} argument, avoiding the need
to regenerate the grid for every estimation.
}
\seealso{
\code{\link{CircularMAD}}
\code{\link{CircularLMS}}
\code{\link{CircularLTS}}
}
\examples{
# Generate a smaller calibration grid
grid <- generate_results_grid_sigma(
  sigma_seq = seq(0.1, 1.5, length.out = 50)
)

head(grid)
}
