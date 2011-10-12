\name{climdex.tn10p}
\alias{climdex.tn10p}
\alias{climdex.tx10p}
\alias{climdex.tn90p}
\alias{climdex.tx90p}
\title{climdex.tn10p}
\description{
  These functions compute the climdex indices TN10p, TX10p, TN90p, and
  TN10p.
}
\usage{
climdex.tn10p(ci)
climdex.tx10p(ci)
climdex.tn90p(ci)
climdex.tx90p(ci)
}
\arguments{
  \item{ci}{Object of type climdexInput.}
}

\details{
  These functions take a climdexInput object as input and compute the
  climdex indices TN10p, TX10p, TN90p, and TN10p. These are:
  \itemize{
    \item{TN10p: }{Monthly percent of values below the 10th percentile
  of baseline daily minimum temperature.}
    \item{TX10p: }{Monthly percent of values below the 10th percentile
  of baseline daily maximum temperature.}
    \item{TN90p: }{Monthly percent of values above the 90th percentile
  of baseline daily minimum temperature.}
    \item{TX90p: }{Monthly percent of values above the 90th percentile
  of baseline daily maximum temperature.}
  }  
  
  Computation of these percentiles involves use of a boostrap procedure,
  described below but described in more depth in [Zhang, 2005]. 
  
  Computation of these values outside of the base period involves
  comparing the temperature data for each day with the corresponding
  percentiles for a 5 day running window surrounding that day. The
  resulting monthly series is then the monthly percentage of values that
  meet the criteria.

  Computation of these values inside the base period is more
  complicated. It involves comparison of the daily temperature data with
  the corresponding day of temperature data in each of (n - 1) sets of
  data. The sets consist of the data for the base period with the
  current year replaced with each of the other years. The results of
  these comparisons are then averaged to give a value between 0 and
  1. Finally, the resulting daily series is aggregated to a monthly
  series by averaging these daily values and multiplying by 100 to give
  a monthly percentile value.
}

\value{
  A vector containing the value of the index for each month.
}

\note{
  These functions may calculate slightly different results than fclimdex.

  The bootstrapping method is not well defined for cases where the base
  data contains numerous missing values.  Because of that, this code
  (and fclimdex) are not very robust against missing values with respect
  to these indicies.  When
  computing percentiles inside the base period, both this implementation
  and fclimdex do not divide through by the number of non-missing values
  when aggregating the values inside the base period. Instead, they
  divide through by the number of base years minus one. This will
  result in a negative bias when missing values are present.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.
## tmax.data, tmin.data, and prec.data are vectors of real numbers.
\dontrun{ci <- climdexInput.raw(tmax.data, tmin.data, prec.data,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))}

## Create a monthly timeseries of percentage of daily minimum
## temperature values which fall below the 10th percentile.
\dontrun{tn10p <- climdex.tn10p(ci)}

## Create a monthly timeseries of percentage of daily maximum
## temperature values which fall below the 10th percentile.
\dontrun{tx10p <- climdex.tx10p(ci)}

## Create a monthly timeseries of percentage of daily minimum
## temperature values which are above the 90th percentile.
\dontrun{tn90p <- climdex.tn90p(ci)}

## Create a monthly timeseries of percentage of daily maximum
## temperature values which are above the 90th percentile.
\dontrun{tx90p <- climdex.tx90p(ci)}
}

\seealso{
  \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}

  Zhang, X., 2005: Avoiding inhomogeneity in percentile-based indices of
  temperature extremes. Journal of Climate 18.11 (2005):1641-.
}  

\keyword{ts}
\keyword{climate}