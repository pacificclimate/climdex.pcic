\name{climdex.su}
\alias{climdex.fd}
\alias{climdex.su}
\alias{climdex.id}
\alias{climdex.tr}
\title{climdex.su}
\description{
  These functions compute the climdex indices FD, SU, ID, and TR.
}
\usage{
climdex.fd(ci)
climdex.su(ci)
climdex.id(ci)
climdex.tr(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices FD, SU, ID, and TR. These are:
  \itemize{
    \item{FD: }{Number of frost days; annual count of days days where
      daily minimum temperature drops below 0 degrees Celsius.}
    \item{SU: }{Number of summer days; annual count of days where daily
      maximum temperature exceeds 25 degrees Celsius.}
    \item{ID: }{Number of icing days; annual count of days where daily
      maximum temperature is below 0 degrees Celsius.}
    \item{TR: }{Number of tropical nights; annual count of days where
      daily minimum temperature exceeds 20 degrees Celsius.}
  }  
}

\note{
  These functions may create slightly different results than fclimdex.
}

\value{
  A vector containing the number of days meeting the criteria of the
  index for each year.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create an annual timeseries of the number of frost days.
\dontrun{fd <- climdex.fd(ci)}

## Create an annual timeseries of the number of summer days.
\dontrun{su <- climdex.su(ci)}

## Create an annual timeseries of the number of icing days.
\dontrun{id <- climdex.id(ci)}

## Create an annual timeseries of the number of tropical nights.
\dontrun{tr <- climdex.tr(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}