\name{climdex.rx5day, climdex.rx1day}
\alias{climdex.rx5day}
\alias{climdex.rx1day}
\title{Monthly Maximum Consecutive 5- and 1-day Precipitation}
\description{
  These functions compute the climdex indices Rx5day and Rx1day.
}
\usage{
climdex.rx5day(ci, freq=c("monthly", "annual"), center.mean.on.last.day=FALSE)
climdex.rx1day(ci, freq=c("monthly", "annual"))
}
\arguments{
  \item{ci}{Object of type climdexInput.}
  \item{freq}{Time frequency to aggregate to.}
  \item{center.mean.on.last.day}{Whether to center the 5-day running
    mean on the last day of the window, instead of the center day.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices Rx5day and Rx1day. These are:
  \itemize{
    \item{Rx5day: }{Monthly or annual maximum 5-day consecutive precipitation.}
    \item{Rx1day: }{Monthly or annual maximum 1-day precipitation.}
  }  
}

\value{
  A vector containing the value of the index for each month of each year.
}

\note{
  The default behaviour of climdex.rx5day differs somewhat from
  fclimdex, as fclimdex and climdex.pcic differ on the definition of
  Rx5day. The running sum series computed by fclimdex is off by 2 days,
  and the first day a running sum can be computed for is left out
  entirely. The behaviour of fclimdex can be replicated by setting
  center.mean.on.last.day to TRUE.
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

## Create a timeseries of monthly maximum 5-day consecutive
## precipitation.
rx5day <- climdex.rx5day(ci)

## Create a timeseries of monthly maximum 1-day precipitation.
rx1day <- climdex.rx1day(ci)

## Print these out for testing purposes.
rx5day
rx1day
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}