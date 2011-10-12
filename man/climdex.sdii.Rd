\name{climdex.sdii}
\alias{climdex.sdii}
\title{climdex.sdii}
\description{
  This function computes the climdex index SDII.
}
\usage{
climdex.sdii(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  \code{climdex.sdii} computes the climdex index SDII, or Simple
  Precipitation Intensity Index. This is defined as the sum of
  precipitation in wet days (days with preciptitation over 1mm) during
  the year divided by the number of wet days in the year.
}

\value{
  A vector containing the value of the index for each year.
}

\note{
  fclimdex rounds to 1 decimal place, whereas climdex.sdii does
  not. This results in some small differences.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create a timeseries of annual SDII values.
\dontrun{sdii <- climdex.sdii(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}