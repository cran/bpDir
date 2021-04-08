\name{CircularBoxplot}
\alias{CircularBoxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Circular Box-and-Wisker Plot
%%  ~~function to do ... ~~
}
\description{
The CircularBoxplot function produces a box-and-whisker-plot for circular data.
%%  ~~ A concise (1-5 lines) description of what the function does.
}
\usage{
CircularBoxplot ( A, template="degrees", place="none", units="degrees", marg=
                   "large", shrink = 1.5,  H=FALSE, stack=FALSE, constant=
                   "optimal")

}
%- maybe also 'usage' for other objects documented here.
\arguments{


    \item{A}{numeric; a vector of class circular in \code{units=c("degrees", "radians")}}

%  \item{method}{character; the type of depth to be used to produce an outer-inward ranking of \code{A}. It can be set to \code{"Tukey"}(default), \code{"simplicial"}, \code{"halfspace", \code{"arcdistance"}}}

\item{template}{character; graphical parameter to  set the template to be used in the plot. Options are \code{"degrees"}(default), \code{"radians"}  and  \code{"geographics"}. If \code{NULL}, a generic template can be defined, with angles displayed according to the graphical parameter  \code{place}.}

%  \item{axes}{logical; graphical parameter. If \code{FALSE} axes labels are not displayed.}

  \item{place}{character; graphical parameter. Active only when \code{template = NULL} to add an additional generic template either \code{"outside"} or \code{"inside"} the circle. \cr If \code{place = "none"}(default), no additional template is plotted.}

\item{units}{character; graphical parameter to select the unit of measure to be displayed on the axes when \code{template = NULL} and \code{place=c("inside","outside")}. Options are \code{"degrees"}(default) or \code{"radians"}.}

\item{marg}{character; graphical parameter to set the plot region to either \code{"large"}(default) or \code{"small"}, expressed in number of lines of text.}

  \item{shrink}{numeric; graphical parameter to control the size of the plotted circle. Default is 1.5. Larger values shrink the circle, while smaller values enlarge the circle.}

  \item{H}{logical; graphical parameter to specify if observations between the hinges and the whiskers are plotted. Default is \code{FALSE}.}

  \item{stack}{If \code{TRUE} and H=\code{FALSE}, then far out values are stacked on the circumference.}

  \item{constant}{numeric;  the boxplot multiplying factor (or "constant") that determines how far the fences extend out from the box. When \code{"optimal"}(default), it ensures the probability that an observation falls outside the fences (i.e., it will be a far out value) is 0,7\% under a von Mises distribution (see Web Appendix A in Buttarazzi et al., 2018).}


%%     ~~Describe \code{x} here~~
}
\details{The arguments \code{stack}, \code{shrink} are parameters passed to \code{plot.circular.}

%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{farout}{a vector of class circular in radians, containing the values of the
points lying outside the whiskers.}

\item{statistics}{a dataframe of summary statistics in degrees: circular median, hinges and whiskers.}

\item{constant}{the value of the constant used to draw the boxplot.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Buttarazzi D., Pandolfo G., Porzio G.C. (2018). A boxplot for circular data, Biometrics.

%% ~put references to the literature/web site here ~
}
\author{Davide Buttarazzi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{plot.circular, points.circular,
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
# Circular Boxplot on Vanishing directions of homing pigeons
#install.packages("circular")
#install.packages("plotrix")
require(circular)
require(plotrix)
#help(fisherB12c)
data(fisherB12c)
CircularBoxplot(fisherB12c, template="geographics")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{device}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
