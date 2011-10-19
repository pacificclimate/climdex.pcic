\name{climdex.r95ptot, climdex.r99ptot}
\alias{climdex.r95ptot}
\alias{climdex.r99ptot}
\alias{climdex.prcptot}
\title{Total Daily Precipitation Exceeding 95\%ile, 99\%ile Threshold}
\description{
  These functions compute the climdex indices R95pTOT, R99pTOT, and
  PRCPTOT.
}
\usage{
climdex.r95ptot(ci)
climdex.r99ptot(ci)
climdex.prcptot(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices R95pTOT, R99pTOT, and PRCPTOT. These are:
  \itemize{
    \item{R95pTOT: }{Annual sum of precipitation in days where daily
  precipitation exceeds the 95th percentile of daily precipitation in
  the base period.}
    \item{R99pTOT: }{Annual sum of precipitation in days where daily
  precipitation exceeds the 99th percentile of daily precipitation in
  the base period.}
    \item{PRCPTOT: }{Annual sum of precipitation in days where daily
  precipitation is at least 1mm/day.}
  }  
}

\value{
  A vector containing an annual timeseries of precipitation exceeding
  the appropriate threshold.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.

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

## Create an annual timeseries of the sum of precipitation where
## precipitation exceeds the 95th percentile of precipitation in the
## base period.
r95ptot <- climdex.r95ptot(ci)

## Create an annual timeseries of the sum of precipitation where
## precipitation exceeds the 99th percentile of precipitation in the
## base period.
r99ptot <- climdex.r99ptot(ci)

## Create an annual timeseries of the sum of precipitation where
## precipitation is at least 1mm/day.
prcptot <- climdex.prcptot(ci)

## Print these out for testing purposes.
r95ptot
r99ptot
prcptot
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}