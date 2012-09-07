\name{nday.consec.prec.max}
\alias{nday.consec.prec.max}
\title{Number of days (less than, greater than, etc) a threshold}
\description{
  Produces sums of values that exceed (or are below) the specified threshold.
}
\usage{
nday.consec.prec.max(daily.prec, date.factor, ndays, center.mean.on.last.day=FALSE)
}
\arguments{
  \item{daily.prec}{Daily timeseries of precipitation.}
  \item{date.factor}{Factor to aggregate by.}
  \item{ndays}{Number of days in the running window.}
  \item{center.mean.on.last.day}{Whether to center the n-day running
    mean on the last day of the series, instead of the middle day.}
}
\details{
  This function takes a data series, the number of days in the running
  window, a date factor to aggregate by, and an optional modifier
  parameter (center.mean.on.last.day). It computes the n-day running sum
  of precipitation and returns the maximum n-day total precipitation per
  unit time, as defined by \code{date.factor}.
}
\value{
  A vector consisting of the maximum n-day sum of precipitation per time
  interval.
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

## Compute rx5day on a monthly basis.
rx5day <- nday.consec.prec.max(ci@prec, ci@monthly.factor, 5)
}
\keyword{ts}
\keyword{climate}