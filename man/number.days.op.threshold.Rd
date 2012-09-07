\name{number.days.op.threshold}
\alias{number.days.op.threshold}
\title{Number of days (less than, greater than, etc) a threshold}
\description{
  Produces sums of values that exceed (or are below) the specified threshold.
}
\usage{
number.days.op.threshold(temp, date.factor, threshold, op="<")
}
\arguments{
  \item{temp}{Sequence temperature values.}
  \item{date.factor}{Factor to aggregate by.}
  \item{threshold}{Threshold to use.}
  \item{op}{Operator to use for comparison.}
}
\details{
  This function takes a data series, a threshold, an operator, and a
  factor to aggregate by. It uses the operator to compare the threshold
  to the data series, creating a series of booleans, then sums the
  booleans according to the factor.
}
\value{
  A vector consisting of the number of values that meet the criteria in
  the given time period (as specified by \code{date.factor}).
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

## Calculate frost days.
fd <- number.days.op.threshold(ci@tmin, ci@annual.factor, 0, "<")
}
\keyword{ts}
\keyword{climate}