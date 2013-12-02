\name{spell.length.max}
\alias{spell.length.max}
\title{Maximum spell length}
\description{
  This function returns the longest string of days which exceed or are
  below the given threshold.
}
\usage{
spell.length.max(daily.prec, date.factor, threshold, op, spells.can.span.years)
}
\arguments{
  \item{daily.prec}{Data to compute index on.}
  \item{date.factor}{Date factor to split by.}
  \item{threshold}{The threshold to compare to.}
  \item{op}{The operator to use to compare data to threshold.}
  \item{spells.can.span.years}{Whether spells can span years.}
}
\details{
  This routine compares data to the threshold using the given operator,
  generating a series of TRUE or FALSE values. It then computes the
  lengths of sequences of TRUE values (spells) and chooses the longest
  spell in each period (as defined by date.factor).

  The \code{spells.can.span.years} option controls whether spells must
  always terminate at the end of a period, or whether they may continue
  until the criteria ceases to be met or the end of the data is
  reached. The default for fclimdex is TRUE.
}
\value{
  A timeseries of maximum spell lengths for each period.
}
\examples{
prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
phony.date.factor <- factor(rep(1:2, each=5))

## With spells spanning years...
cwd <- spell.length.max(prec.dat, phony.date.factor, 1, ">=", TRUE)

## Without spells spanning years...
altcwd <- spell.length.max(prec.dat, phony.date.factor, 1, ">=", FALSE)
}
\seealso{
  \code{\link{climdex.cdd}}.
}
\keyword{ts}
\keyword{climate}