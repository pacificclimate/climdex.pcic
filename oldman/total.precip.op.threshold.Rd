\name{total.precip.op.threshold}
\alias{total.precip.op.threshold}
\title{Sum of precipitation above a threshold}
\description{
  This function returns the sum of values above a threshold for each
  period (as defined by date.factor).
}
\usage{
total.precip.op.threshold(daily.prec, date.factor, threshold, op)
}
\arguments{
  \item{daily.prec}{Data to compute index on.}
  \item{date.factor}{Date factor to split by.}
  \item{threshold}{The threshold to compare to.}
  \item{op}{The operator to use to compare data to threshold.}
}
\details{
  This routine sums up all values which exceed or are below (depending
  on op) the given threshold.
}
\value{
  A timeseries of sums of numbers above the threshold for each period.
}
\examples{
prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
phony.date.factor <- factor(rep(1:2, each=5))

## Compute equiv of PRCPTOT
prec.sum <- total.precip.op.threshold(prec.dat, phony.date.factor, 1, ">=")
}
\seealso{
  \code{\link{climdex.r99ptot}}.
}
\keyword{ts}
\keyword{climate}