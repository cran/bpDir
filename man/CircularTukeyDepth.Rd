\name{CircularTukeyDepth}
\alias{CircularTukeyDepth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Circular Tukey Depth
}
\description{
Empirical depth ranking based on Tukey’s ranking procedure (Tukey, 1977) adapted to the circular setting (Buttarazzi et al, 2018).
}
\usage{
CircularTukeyDepth(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric; a vector of class \code{circular} in \code{units=c("degrees", "radians")}}
}

\value{
\item{depth}{a two column matrix of circular observations with corresponding depth value}
\item{median}{the depth-based circular median as an object of class \code{circular}}
\item{iqr}{the depth-based circular inter-quartile range as an object of class \code{circular}}
}

\references{

Buttarazzi D., Pandolfo G., Porzio G.C. (2018). A boxplot for circular data, Biometrics.

Tukey, J. W. (1977). Exploratory Data Analysis. MA: AddisonWesley, Reading (UK).

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
CircularTukeyDepth(data1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ depth }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
