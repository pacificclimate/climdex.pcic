\name{threshold.exceedance.duration.index}
\alias{threshold.exceedance.duration.index}
\title{Sum of spell lengths exceeding daily threshold}
\description{
  This function returns the number of spells of more than
  \code{min.length} days which exceed or are below the given threshold.
}
\usage{
threshold.exceedance.duration.index(daily.temp, date.factor, jdays, thresholds, op=">", min.length=6, spells.can.span.years=TRUE)
}
\arguments{
  \item{daily.temp}{Data to compute index on.}
  \item{date.factor}{Date factor to split by.}
  \item{jdays}{Timeseries of days of year.}
  \item{thresholds}{The thresholds to compare to.}
  \item{op}{The operator to use to compare data to threshold.}
  \item{min.length}{The minimum spell length to be considered.}
  \item{spells.can.span.years}{Whether spells can span years.}
}
\details{
  This routine compares data to the thresholds using the given operator,
  generating a series of TRUE or FALSE values; these values are then
  filtered to remove any sequences of less than \code{min.length} days
  of TRUE values. It then computes the lengths of the remaining
  sequences of TRUE values (spells) and sums their lengths.

  The \code{spells.can.span.years} option controls whether spells must
  always terminate at the end of a period, or whether they may continue
  until the criteria ceases to be met or the end of the data is
  reached. The default for fclimdex is FALSE.
}
\value{
  A timeseries of maximum spell lengths for each period.
}
\examples{
prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
phony.date.factor <- factor(rep(1:2, each=5))

## With spells spanning years...
alttedi <- threshold.exceedance.duration.index(prec.dat,
phony.date.factor, rep(1:5, 2), rep(1, 5), ">=", 2, TRUE)

## Without spells spanning years...
tedi <- threshold.exceedance.duration.index(prec.dat, phony.date.factor,
rep(1:5, 2), rep(1, 5), ">=", 2, FALSE)
}
\seealso{
  \code{\link{climdex.wsdi}}.
}
\keyword{ts}
\keyword{climate}