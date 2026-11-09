\name{Spring}
\alias{Spring}
\title{Spring wind directions}
\description{
Wind-direction measurements for spring, stored as circular observations.
}
\usage{
data("Wind_directions")
}
\format{
A matrix of class \code{circular}, with angles in radians. The supplied object includes a column label \code{Spring2010}.
}
\details{
This dataset was supplied in the file \file{Wind_directions.RData}.
The observations use the \code{circular} class from the \pkg{circular} package.
}
\source{
Demni et al. (2026), Notched violin plots for circular data. In: Martella, F., Arima, S., Marino, M.F., Mollica, C. (eds) Statistical Science: From Theory to Applied Research I. SIS-FENStatS 2026. Italian Statistical Society Series on Advances in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-032-30873-3 10.
}
\seealso{
\code{\link[circular]{circular}}
}
\examples{
data("Wind_directions")
str(Spring)
summary(Spring)
}
\keyword{datasets}
