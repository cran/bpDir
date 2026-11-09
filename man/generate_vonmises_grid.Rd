\name{generate_vonmises_grid}
\alias{generate_vonmises_grid}
\title{Generate a Calibration Grid for the von Mises Concentration Parameter}
\description{
Generates a calibration grid relating robust measures of circular
dispersion to the circular standard deviation (CSD) and the
concentration parameter (\eqn{\kappa}) of a von Mises distribution.

The grid is used by the estimation procedure to convert an observed
robust circular dispersion statistic into an estimate of CSD through
monotone spline interpolation, and subsequently into an estimate of
the von Mises concentration parameter.
}
\usage{
generate_vonmises_grid(
  rho_seq = seq(0.01, 0.99, by = 0.002)
)
}
\arguments{
  \item{rho_seq}{
    A numeric vector of values of the mean resultant length
    \eqn{\rho} for which the calibration grid is generated.
    Values must lie strictly between 0 and 1.
    The default sequence ranges from 0.01 to 0.99 in increments
    of 0.002.
  }
}
\value{
A data frame with the following columns:
\describe{
  \item{kappa}{
    Values of the von Mises concentration parameter, calculated
    by applying \code{circular::A1inv} to the corresponding values
    of \eqn{\rho}.
  }
  \item{rho}{
    Values of the mean resultant length.
  }
  \item{CSD}{
    The circular standard deviation, calculated as
    \eqn{2(1-\rho)}.
  }
  \item{CMAD}{
    Calibration values obtained using the 0.75 quantile of the
    centred von Mises distribution, for the CMAD-based estimation
    procedure.
  }
  \item{CLMS}{
    Calibration values obtained using the 0.75 quantile of the
    centred von Mises distribution, for the CLMS-based estimation
    procedure.
  }
  \item{CLTS}{
    Calibration values computed using the internal procedure
    for the von Mises distribution, for the CLTS-based estimation
    procedure.
  }
}
}
\details{
The calibration grid is constructed by evaluating the dispersion
statistics over a sequence of von Mises distributions with different
values of the concentration parameter \eqn{\kappa}.

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

# Generate a custom calibration grid
custom_grid <- generate_vonmises_grid(
  rho_seq = seq(0.05, 0.95, length.out = 50)
)

head(custom_grid)
}

