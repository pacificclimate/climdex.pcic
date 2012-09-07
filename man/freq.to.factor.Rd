\name{freq.to.factor}
\alias{freq.to.factor}
\alias{freq.to.namask}
\title{Get NA masks or date factors given frequency definition}
\description{
  These functions take a frequency definition and return data
  appropriate for that output frequency.
}
\usage{
freq.to.factor(ci, freq)
freq.to.namask(ci, freq)
}
\arguments{
  \item{ci}{ClimdexInput object.}
  \item{freq}{Frequency (either monthly or annual).}
}
\details{
  These are convenience functions which simply return the appropriate NA
  mask or date factor given an output time frequency (monthly or annual).
}
\value{
  An NA mask or date factor for the given frequency.
}
\examples{
## Parse the dates into PCICt.
tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
"jday")]), format="\%Y \%j", cal="gregorian")
tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
"jday")]), format="\%Y \%j", cal="gregorian")
prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
"jday")]), format="\%Y \%j", cal="gregorian")

## Load the data in.
ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))

## Get appropriate NA mask and date factor
na.mask.ann <- freq.to.namask(ci, "annual")
date.factor.ann <- freq.to.factor(ci, "annual")
}
\keyword{ts}
\keyword{climate}