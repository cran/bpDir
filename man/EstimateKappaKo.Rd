\name{EstimateKappaKo}
\alias{EstimateKappaKo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust estimation of the concentration parameter of the von Mises-Fisher distribution
}
\description{
Robust estimation of the concentration parameter of the von Mises-Fisher distribution (Ko, 1992).
}
\usage{
EstimateKappaKo(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{numeric; a vector of class \code{circular} in \code{units=c("radians")}}
}

\value{
A numeric value
}

\references{
Ko, D. (1992), Robust estimation of the concentration parameter of the von Mises-Fisher distribution.
Ann. Stat. 20(2) 917-928

%% ~put references to the literature/web site here ~
}
\author{Davide Buttarazzi
%%  ~~who you are~~
}

\seealso{CircularBoxplot}

\examples{
#install.packages("circular")
require(circular)
data1 <- rvonmises(50, circular(pi/2), 5)
EstimateKappaKo(data1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circular }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
