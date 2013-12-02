\name{select.blocks.gt.length}
\alias{select.blocks.gt.length}
\title{Select blocks of TRUE values of sufficient length.}
\description{
  Produces a sequence of booleans of the same length as input, with
  sequences of TRUE values shorter than n replaced with FALSE.
}
\usage{
select.blocks.gt.length(d, n, na.value=FALSE)
}
\arguments{
  \item{d}{Sequence of booleans.}
  \item{n}{Longest sequence of TRUE to replace with FALSE.}
  \item{na.value}{Values to replace NAs with.}
}
\details{
  This function takes a series of booleans and returns a sequence of
  booleans of equal length, with all sequences of TRUE of length
  \code{n} or shorter replaced with sequences of FALSE. NA values are
  replaced with \code{na.value}.
}
\value{
  A vector of booleans, with the length \code{n} or less sequences of
  TRUE replaced with FALSE.
}
\examples{
## Return only the first sequence of TRUE... second sequence will be FALSE.
foo <- select.blocks.gt.length(c(rep(TRUE, 4), FALSE, rep(TRUE, 3)), 3)
}
\keyword{ts}
\keyword{climate}