\name{climdex.r95ptot}
\alias{climdex.r95ptot}
\alias{climdex.r99ptot}
\alias{climdex.prcptot}
\title{climdex.r95ptot}
\description{
  These functions compute the climdex indices R95pTOT, R99pTOT, and
  PRCPTOT.
}
\usage{
climdex.r95ptot(ci)
climdex.r99ptot(ci)
climdex.prcptot(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices R95pTOT, R99pTOT, and PRCPTOT. These are:
  \itemize{
    \item{R95pTOT: }{Annual sum of precipitation in days where daily
  precipitation exceeds the 95th percentile of daily precipitation in
  the base period.}
    \item{R99pTOT: }{Annual sum of precipitation in days where daily
  precipitation exceeds the 99th percentile of daily precipitation in
  the base period.}
    \item{PRCPTOT: }{Annual sum of precipitation in days where daily
  precipitation is at least 1mm/day.}
  }  
}

\value{
  A vector containing an annual timeseries of precipitation exceeding
  the appropriate threshold.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create an annual timeseries of the sum of precipitation where
## precipitation exceeds the 95th percentile of precipitation in the
## base period.
\dontrun{r95ptot <- climdex.r95ptot(ci)}

## Create an annual timeseries of the sum of precipitation where
## precipitation exceeds the 99th percentile of precipitation in the
## base period.
\dontrun{r99ptot <- climdex.r99ptot(ci)}

## Create an annual timeseries of the sum of precipitation where
## precipitation is at least 1mm/day.
\dontrun{prcptot <- climdex.prcptot(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}