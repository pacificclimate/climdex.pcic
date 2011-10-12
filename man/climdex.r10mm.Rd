\name{climdex.r10mm}
\alias{climdex.r10mm}
\alias{climdex.r20mm}
\alias{climdex.rnnmm}
\title{climdex.r10mm}
\description{
  These functions compute the climdex indices R10mm, R20mm, and Rnnmm.
}
\usage{
climdex.r10mm(ci)
climdex.r20mm(ci)
climdex.rnnmm(ci, threshold)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
  \item{threshold}{The threshold to be used for Rnnmm.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices R10mm, R20mm, and Rnnmm. These are:
  \itemize{
    \item{R10mm: }{Annual count of days where daily precipitation is more
  than 10mm.}
    \item{R20mm: }{Annual count of days where daily precipitation is more
  than 20mm.}
    \item{Rnnmm: }{Annual count of days where daily precipitation is more
  than the number of millimeters specified by \code{threshold}.}
  }  
}

\value{
  A vector containing the value of the index for each year.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create an annual timeseries of the number of days where
## precipitation exceeds 10mm/day.
\dontrun{r10mm <- climdex.r10mm(ci)}

## Create an annual timeseries of the number of days where
## precipitation exceeds 20mm/day.
\dontrun{r20mm <- climdex.r20mm(ci)}

## Create an annual timeseries of the number of days where
## precipitation exceeds 25mm/day.
\dontrun{r25mm <- climdex.rnnmm(ci, 25)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}