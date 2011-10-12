\name{climdex.rx5day}
\alias{climdex.rx5day}
\alias{climdex.rx1day}
\title{climdex.rx5day}
\description{
  These functions compute the climdex indices Rx5day and Rx1day.
}
\usage{
climdex.rx5day(ci)
climdex.rx1day(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices Rx5day and Rx1day. These are:
  \itemize{
    \item{Rx5day: }{Monthly maximum 5-day consecutive precipitation.}
    \item{Rx1day: }{Monthly maximum 1-day precipitation.}
  }  
}

\value{
  A vector containing the value of the index for each month of each year.
}

\note{
  Behaviour of climdex.rx5day differs somewhat from fclimdex, as
  fclimdex contains a bug. The running sum series computed by fclimdex
  is off by 2 days, and the first day a running sum can be computed for
  is left out entirely. This can result in wet periods near a month
  boundary going into a different month.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create a timeseries of monthly maximum 5-day consecutive
## precipitation.
\dontrun{rx5day <- climdex.rx5day(ci)}

## Create a timeseries of monthly maximum 1-day precipitation.
\dontrun{rx1day <- climdex.rx1day(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}