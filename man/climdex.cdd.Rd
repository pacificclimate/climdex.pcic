\name{climdex.cdd, climdex.cwd}
\alias{climdex.cdd}
\alias{climdex.cwd}
\title{Maximum Consecutive Wet and Dry Days}
\description{
  These functions compute the climdex indices CDD and CWD.
}
\usage{
climdex.cdd(ci)
climdex.cwd(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices CDD and CWD. These are:
  \itemize{
    \item{CDD: }{Annual maximum length of dry spell, in days.}
    \item{CWD: }{Annual maximum length of wet spell, in days.}
  }  

  Wet spells are considered to be sequences of days where daily
  precipitation is at least 1mm per day.

  Dry spells are considered to be sequences of days where daily
  preciptation is less than 1mm per day.
}

\value{
  A vector containing the length of the dry or wet spell for each year.
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

## Create an annual timeseries of the number of consecutive days where
## precipitation was less than 1mm/day.
cdd <- climdex.cdd(ci)

## Create an annual timeseries of the number of consecutive days where
## precipitation was at least 1mm/day.
cwd <- climdex.cwd(ci)

## Print these out for testing purposes.
cdd
cwd
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}