\name{climdex.gsl}
\alias{climdex.gsl}
\title{Growing Season Length}
\description{
  This function computes the growing season length (GSL) given the
  input.
}
\usage{
climdex.gsl(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  This function takes a climdexInput object as input and computes the
  growing season length based on this data.

  Growing season length as defined by the climdex indices is the number
  of days between the start of the first spell of warm days in the first
  half of the year, and the start of the first spell of cold days in the
  second half of the year. Spells of warm days are defined as six or
  more days with mean temperature above 5 degrees Celsius; spells of
  cold days are defined as six or more days with a mean temperature
  below 5 degrees Celsius.
}

\note{
  \code{climdex.gsl} only implements GSL for the northern
  hemisphere. Should you want to use this function in the southern
  hemisphere, we suggest you move the data forward or backwards by six
  months.
  
  Note that fclimdex contains a number of bugs with growing season
  length. These include:
  \itemize{
    \item{Reporting 1 less day of growing season than it should.}
    \item{Reporting a maximum of n-6 days of growing season in the
      tropics.}
    \item{Invalid results in the southern hemisphere overall.}
  }

  Because of these issues, do not expect identical results from fclimdex
  and climdex.pcic .
}

\value{
  A vector containing the number of days meeting the criteria of the
  index for each year.
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
gsl <- climdex.gsl(ci)

## Print these out for testing purposes.
gsl
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}