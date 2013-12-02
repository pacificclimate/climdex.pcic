\name{climdex.quantile}
\alias{climdex.quantile}
\title{Climdex quantile function}
\description{
  This function implements R's type=8 in a more efficient manner.
}
\usage{
climdex.quantile(x, q=c(0, 0.25, 0.5, 0.75, 1))
}
\arguments{
  \item{x}{Data to compute quantiles on.}
  \item{q}{Quantiles to be computed.}
}
\details{
  This is a reimplementation of R's type=8 created to improve the
  efficiency of this package.
}
\value{
  A vector of the quantiles in question.
}
\examples{
## Compute 10th, 50th, and 90th percentile of example data.
climdex.quantile(1:10, c(0.1, 0.5, 0.9))
}
\seealso{
  \code{\link{quantile}}
}
\keyword{ts}
\keyword{climate}