\name{climdex.txx, climdex.tnx, climdex.txn, climdex.tnn}
\alias{climdex.txx}
\alias{climdex.tnx}
\alias{climdex.txn}
\alias{climdex.tnn}
\title{Monthly Maximum/Minimum of Daily Maximum/Minimum Temperature}
\description{
  These functions compute the climdex indices TXx, TNx, TXn, and TNn.
}
\usage{
climdex.txx(ci)
climdex.tnx(ci)
climdex.txn(ci)
climdex.tnn(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices TXx, TNx, TXn, and TNn. These are:
  \itemize{
    \item{TXx: }{Monthly maximum of daily maximum temperature.}
    \item{TNx: }{Monthly maximum of daily minimum temperature.}
    \item{TXn: }{Monthly minimum of daily maximum temperature.}
    \item{TNn: }{Monthly minimum of daily minimum temperature.}
  }  
}

\value{
  A vector containing the value of the index for each month.
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

## Create a monthly timeseries of maximum daily maximum temperature.
txx <- climdex.txx(ci)

## Create a monthly timeseries of maximum daily minimum temperature.
tnx <- climdex.tnx(ci)

## Create a monthly timeseries of minimum daily maximum temperature.
txn <- climdex.txn(ci)

## Create a monthly timeseries of minimum daily minimum temperature.
tnn <- climdex.tnn(ci)

## Print these out for testing purposes.
txx
tnx
txn
tnn
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}