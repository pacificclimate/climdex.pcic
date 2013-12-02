\name{get.series.lengths.at.ends}
\alias{get.series.lengths.at.ends}
\title{Get series length at ends}
\description{
  This function takes a series of boolean values and returns a list of
  integers of the same length corresponding to the lengths at the ends
  of sequences of TRUE values.
}
\usage{
get.series.lengths.at.ends(x, na.value=FALSE)
}
\arguments{
  \item{x}{Sequence of booleans.}
  \item{na.value}{Value to replace NAs with.}
}
\details{
  It can often be useful to know how long a series of boolean values
  is. This function provides a method of knowing where and how long such
  sequences are.
}
\value{
  A vector consisting of the lengths of sequences of TRUE values at the
  location of the last TRUE value in the sequence, and zeroes elsewhere.
}
\examples{
## Get lengths of sequences of TRUE values in a sequence
series.lengths <- get.series.lengths.at.ends(c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE,
TRUE, TRUE, FALSE))

}
\keyword{ts}
\keyword{climate}