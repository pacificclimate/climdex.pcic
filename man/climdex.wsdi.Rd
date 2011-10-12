\name{climdex.wsdi}
\alias{climdex.wsdi}
\alias{climdex.csdi}
\title{climdex.wsdi}
\description{
  These functions compute the climdex indices WSDI and CSDI.
}
\usage{
climdex.wsdi(ci)
climdex.csdi(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices WSDI and CSDI. These are:
  \itemize{
    \item{WSDI: }{Warm Spell Duration Index.}
    \item{CSDI: }{Cold Spell Duration Index.}
  }  

  The warm spell duration index is defined as the number of days each
  year which are part of a "warm spell". A "warm spell" is defined as a
  sequence of 6 or more days in which the daily maximum temperature
  exceeds the 90th percentile of daily maximum temperature for a 5-day
  running window surrounding this day during the baseline period.
  
  Likewise, the cold spell duration index is defined as the number of
  days each year which are part of a "cold spell". A "cold spell" is
  defined as a sequence of 6 or more days in which the daily minimum
  temperature is below the 10th percentile of daily minimum temperature
  for a 5-day running window surrounding this day during the baseline
  period.
}

\value{
  A vector containing the value of the index for each year.
}

\note{
  These functions may calculate slightly different results than fclimdex.

  Behaviour of climdex.wsdi and climdex.csdi differ somewhat from
  fclimdex. fclimdex considers all days in
  a warm or cold spell to be part of the year in which the spell ended.
  climdex.wsdi and climdex.csdi split the spell such that days in each
  spell are allocated to the separate years in the days occurred.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create an annual timeseries of the warm spell duration index.
\dontrun{wsdi <- climdex.wsdi(ci)}

## Create an annual timeseries of the cold spell duration index.
\dontrun{csdi <- climdex.csdi(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}