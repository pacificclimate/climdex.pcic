\name{climdex.txx}
\alias{climdex.txx}
\alias{climdex.tnx}
\alias{climdex.txn}
\alias{climdex.tnn}
\title{climdex.txx}
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
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create a monthly timeseries of maximum daily maximum temperature.
\dontrun{txx <- climdex.txx(ci)}

## Create a monthly timeseries of maximum daily minimum temperature.
\dontrun{tnx <- climdex.tnx(ci)}

## Create a monthly timeseries of minimum daily maximum temperature.
\dontrun{txn <- climdex.txn(ci)}

## Create a monthly timeseries of minimum daily minimum temperature.
\dontrun{tnn <- climdex.tnn(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}