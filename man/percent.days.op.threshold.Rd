\name{percent.days.op.threshold}
\alias{percent.days.op.threshold}
\title{Lengths of strings of TRUE values}
\description{
  Computes fraction of days above or below the baseline threshold for
  each day, and averages them using the date factor passed in.
}
\usage{
percent.days.op.threshold(temp, dates, jdays, date.factor, threshold.outside.base, base.thresholds, base.range, op='<')
}
\arguments{
  \item{temp}{Sequence of temperature values.}
  \item{dates}{Sequence of associated dates.}
  \item{jdays}{Sequence of associated days of year.}
  \item{date.factor}{Factor to aggregate data using.}
  \item{threshold.outside.base}{Sequence of thresholds to be used for
    data outside the base period.}
  \item{base.thresholds}{Data structure containing sets of thresholds to
    be used inside the base period; see \link{climdexInput-class}.}
  \item{base.range}{Date range (type PCICt) of the baseline period.}
  \item{op}{Comparison operator to use.}
}
\details{
  This function computes fractions of days above or below baseline
  thresholds for each day, then aggregates them using
  \code{date.factor}. It is used to implement TN/TX 10/90p.
}
\value{
  A vector consisting of the mean fraction of days above or below the
  supplied set of thresholds.
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

## Compute monthly tx90p.
tx90p <- percent.days.op.threshold(ci@tmax, ci@dates, ci@jdays, ci@monthly.factor, ci@running.pctile.notbase$tx90thresh, ci@running.pctile.base$tx90thresh, ci@base.range, ">") * ci@namask.mon$tmax
}
\seealso{
  \link{climdexInput-class}.
}
\keyword{ts}
\keyword{climate}