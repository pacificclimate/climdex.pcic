\name{climdex.gsl}
\alias{climdex.gsl}
\title{Growing Season Length}
\description{
  This function computes the growing season length (GSL) given the
  input.
}
\usage{
climdex.gsl(ci, gsl.mode=c("GSL", "GSL_first", "GSL_max", "GSL_sum"))
}
\arguments{
  \item{ci}{Object of type climdexInput.}
  \item{gsl.mode}{Growing season length method to use.}
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

  The three alternate modes provided ('GSL_first', 'GSL_max', and
  'GSL_sum') are for testing purposes only. They differ considerably
  from the first ('GSL') mode. All of them use a list of growing seasons
  -- here defined as six or more consecutive days with a mean temperature
  greater than or equal to 5 degrees Celsius, followed by either the end
  of the year or six or more consecutive days with a mean temperature
  less than 5 degrees Celsius. 'GSL_first' returns the first growing
  season found; 'GSL_max' returns the longest growing season found; and
  'GSL_sum' returns the total length of all growing seasons found.
}

\note{
  Note that fclimdex results may differ from results using the first
  ('GSL') mode due to bugs in fclimdex. Please ensure you are using the
  latest version of fclimdex, as there have been numerous bug fixes and
  the results should, at this point, match.

  Please do not use the 'GSL_first', 'GSL_max', or 'GSL_sum' modes for
  anything other than testing purposes at this time, nor should you rely
  on this parameter being present in future versions of climdex.pcic.
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
gsl <- climdex.gsl(ci)

## Print these out for testing purposes.
gsl
}

\seealso{
  \code{\link{growing.season.length}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}