\name{growing.season.length}
\alias{growing.season.length}
\title{Flexible GSL function}
\description{
  This function computes the growing season length (GSL) given the
  input, which is allowed to vary considerably from the ETCCDI
  definitions.
}
\usage{
growing.season.length(daily.mean.temp, date.factor, dates,
northern.hemisphere, min.length=6, t.thresh=5,
gsl.mode=c("GSL", "GSL_first", "GSL_max", "GSL_sum"))
}
\arguments{
  \item{daily.mean.temp}{Timeseries of daily mean temperature (in
    degrees C), padded out to end on a year boundary (ie: starts on
    January 1st of some year, ends on December 31st).}
  \item{date.factor}{Factor of the same length as daily.mean.temp that
    divides the timeseries up into years of data.}
  \item{dates}{The corresponding series of dates.}
  \item{northern.hemisphere}{Whether the data is from the northern hemisphere.}
  \item{min.length}{The minimum number of days above or below the
    threshold temperature that defines the start or end of a growing
    season.}
  \item{t.thresh}{The temperature threshold for being considered part of
    a growing season (in degrees C).}
  \item{gsl.mode}{The growing season length mode (ETCCDI mode is "GSL").}
}
\details{
  This function is the function used to implement
  \code{\link{climdex.gsl}}. It's designed to be flexible to allow for
  experimentation and testing of new thresholds and methods.

  If you need to use this code for experimentation in the southern
  hemisphere, you'll need to rip off the climdex.gsl code to rotate the
  year around so that July 1st is treated as January 1st.
  
  See \code{\link{climdex.gsl}} for more information on what \code{gsl.mode} does.
}
\value{
  A vector containing the number of days in the growing season for each year.
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

## Create an annual timeseries of the growing season length in days.
gsl <- growing.season.length(ci@tavg, ci@annual.factor, ci@dates,
                             ci@northern.hemisphere, gsl.mode="GSL") * ci@namask.ann$tavg

## Print these out for testing purposes.
gsl
}

\seealso{
  \code{\link{climdex.gsl}}, \code{\link{climdexInput.csv}}.
}

\keyword{ts}
\keyword{climate}