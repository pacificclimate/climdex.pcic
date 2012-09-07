\name{simple.precipitation.intensity.index}
\alias{simple.precipitation.intensity.index}
\title{Simple Precipitation Intensity Index}
\description{
  This function implements the ETCCDI Simple Precipitation Intensity Index.
}
\usage{
simple.precipitation.intensity.index(daily.prec, date.factor)
}
\arguments{
  \item{daily.prec}{Data to compute index on.}
  \item{date.factor}{Date factor to split by.}
}
\details{
  The simple precipitation intensity index is computed by taking the sum
  of precipitation in wet days (days with >1mm of precipitation), and
  dividing that by the number of wet days in the period. This gives the
  mean precipitation in wet days.
}
\value{
  The mean precipitation in wet days for each period (as defined by date.factor).
}
\examples{
prec.dat <- c(0.1, 3.0, 4.3, 0.9, 1.3, 6.0, 0, 0, 4.0, 1)
phony.date.factor <- factor(rep(1:2, each=5))
sdii <- simple.precipitation.intensity.index(prec.dat, phony.date.factor)
}
\keyword{ts}
\keyword{climate}