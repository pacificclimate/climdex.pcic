\name{climdex.dtr}
\alias{climdex.dtr}
\title{climdex.dtr}
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
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create a monthly timeseries of mean diurnal temperature range.
\dontrun{dtr <- climdex.dtr(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}