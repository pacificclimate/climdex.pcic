\name{climdex.dtr}
\alias{climdex.dtr}
\title{Mean Diurnal Temperature Range}
\description{
  This function computes the diurnal temperature range on a monthly basis.
}
\usage{
climdex.dtr(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  \code{climdex.dtr} computes the mean daily diurnal temperature range
  on a monthly basis.
}

\value{
  A vector containing the mean monthly diurnal temperature range for
  each month.
}

\note{
  This function creates results which may differ in the 3rd decimal
  place from the results from fclimdex.
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

## Create a monthly timeseries of mean diurnal temperature range.
dtr <- climdex.dtr(ci)

## Print these out for testing purposes.
dtr
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}