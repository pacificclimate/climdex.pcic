\name{climdex.cdd}
\alias{climdex.cdd}
\alias{climdex.cwd}
\title{climdex.cdd}
\description{
  These functions compute the climdex indices CDD and CWD.
}
\usage{
climdex.cdd(ci)
climdex.cwd(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices CDD and CWD. These are:
  \itemize{
    \item{CDD: }{Annual maximum length of dry spell, in days.}
    \item{CWD: }{Annual maximum length of wet spell, in days.}
  }  

  Wet spells are considered to be sequences of days where daily
  precipitation is at least 1mm per day.

  Dry spells are considered to be sequences of days where daily
  preciptation is less than 1mm per day.
}

\value{
  A vector containing the length of the dry or wet spell for each year.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create an annual timeseries of the number of consecutive days where
## precipitation was less than 1mm/day.
\dontrun{cdd <- climdex.cdd(ci)}

## Create an annual timeseries of the number of consecutive days where
## precipitation was at least 1mm/day.
\dontrun{cwd <- climdex.cwd(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}