\name{e_parkinson}
\alias{e_parkinson}
\title{Elbow angles for a Parkinson patient}
\description{
Angular gait measurements for a Parkinson parkinson, recorded as left-elbow angles.
}
\usage{
data("GaitData")
}
\format{
A matrix of class \code{circular}, with angles in radians. The supplied object includes a column label \code{left_elbow_angle_parkinson}.
}
\details{
This dataset was supplied in the file \file{GaitData.RData}.
The observations use the \code{circular} class from the \pkg{circular} package.
}
\source{

Pace et al. 2026, Markerless
gait analysis and circular statistics: a case study in Parkinson's disease classification. In: Sintorn, I.M., Pal, U., Lu, S. (eds), Pattern Recognition.
ICPR 2026 International Workshops. Lecture Notes in Computer Science,
vol. 17110. Springer, Cham. DOI: 10.1007/978-3-032-39379-1 4.

Kour, N., Arora, S.: Vision-based gait dataset for gait analysis. In: 2022 2nd Inter-national Conference on Advance Computing andInnovativeTechnologiesinEngi-neering,pp.2163-2168.IEEE(2022)
}
\seealso{
\code{\link[circular]{circular}}
}
\examples{
data("GaitData")
str(e_parkinson)
summary(e_parkinson)
}
\keyword{datasets}
