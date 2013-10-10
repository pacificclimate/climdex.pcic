\name{climdexInput.raw, climdexInput.csv}
\alias{climdexInput.csv}
\alias{climdexInput.raw}
\alias{get.outofbase.quantiles}
\title{Methods for Creating climdexInput Objects}
\description{
  These functions create climdexInput objects for use with the climdex methods.
}
\usage{
climdexInput.raw(tmax, tmin, prec, tmax.dates, tmin.dates,prec.dates,
base.range=c(1961, 1990), n=5, northern.hemisphere=TRUE,
pad.data.with.first.last.values=FALSE, tavg=NULL, tavg.dates=NULL,
temp.quantiles.notbase=NULL, prec.quantiles=NULL, temp.quantiles.base=NULL)

climdexInput.csv(tmax.file, tmin.file, prec.file,
data.columns=list(tmin="tmin", tmax="tmax", prec="prec"),
base.range=c(1961, 1990), na.strings=NULL, cal="gregorian",
date.types=NULL, n=5, northern.hemisphere=TRUE,
pad.data.with.first.last.values=FALSE, tavg.file=NULL,
temp.quantiles.notbase=NULL, prec.quantiles=NULL, temp.quantiles.base=NULL)

get.outofbase.quantiles(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates,
base.range=c(1961, 1990), n=5, pad.data.with.first.last.values=FALSE)
}
\arguments{
  \item{tmax}{Daily maximum temperature data.}
  \item{tmin}{Daily minimum temperature data.}
  \item{prec}{Daily total precipitation data.}
  \item{tavg}{Daily mean temperature data.}
  \item{tmax.dates}{Dates for the daily maximum temperature data.}
  \item{tmin.dates}{Dates for the daily minimum temperature data.}
  \item{prec.dates}{Dates for the daily total precipitation data.}
  \item{tavg.dates}{Dates for the daily mean temperature data.}
  \item{tmax.file}{Name of file containing daily maximum temperature data.}
  \item{tmin.file}{Name of file containing daily minimum temperature data.}
  \item{prec.file}{Name of file containing daily total precipitation data.}
  \item{tavg.file}{Name of file containing daily mean temperature data.}
  \item{data.columns}{Column names for tmin, tmax, and prec data.}
  \item{date.types}{Column names for tmin, tmax, and prec data (see notes).}
  \item{na.strings}{Strings used for NA values; passed to \code{\link{read.csv}}.}
  \item{cal}{The calendar type used in the input files.}
  \item{base.range}{Years to use for the baseline.}
  \item{n}{Number of days to use as window for daily quantiles.}
  \item{northern.hemisphere}{Whether this point is in the northern hemisphere.}
  \item{pad.data.with.first.last.values}{Base data padding option (see notes).}
  \item{temp.quantiles.notbase}{Out of base quantiles (see notes).}
  \item{prec.quantiles}{Precipitation quantiles (see notes).}
  \item{temp.quantiles.base}{In-base quantiles (see notes).}
}

\details{
  These functions take input climate data at daily resolution, and
  produce as output a ClimdexInput data structure. This data structure
  is then passed to each of the routines used to compute the Climdex
  indices. The indices themselves are specified on the webpage cited in
  the references section. \code{climdexInput.csv} should be used when
  your data is already in a CSV file. \code{climdexInput.raw} should be
  used in all other circumstances, as it is much more flexible.
  
  \code{climdexInput.csv} takes as arguments three input filenames: one for
  each of daily maximum temperature, minimum temperature, and total
  preciptation.
  
  \code{climdexInput.raw} takes as arguments the input daily maximum
  temperature, minimum temperature, and total precipitation data along
  with the associated series of dates.

  Daily mean temperature data can optionally be supplied; if it is not
  supplied, it will be computed by taking the mean of the minimum and
  maximum daily temperature for each day.
}

\value{
  An object of class \code{\link{climdexInput-class}} for use with other climdex methods.
}

\note{
  Units are assumed to be mm/day for precipitation and degrees Celsius
  for temperature. No units conversion is performed internally.

  The \code{tmax.dates}, \code{tmin.dates}, and \code{prec.dates}
  arguments are vectors of type \code{PCICt}.
  
  The \code{data.columns} argument is a vector consisting of named items
  tmax, tmin, and prec. These named items are used as the column names
  in their respective files when loading in CSV.
  
  The \code{cal} argument is a textual description of the calendar type,
  as described in the documentation for \code{\link{as.PCICt}}.
  
  The \code{date.types} argument is a list of lists containing two named
  items: \code{fields}, and \code{format}. The \code{fields} item is a
  vector of names consisting of the columns to be concatenated together
  with spaces. The \code{format} item is a date format as taken by
  \code{strptime}.

  The \code{base.range} argument is a pair of 4 digit years which bound
  the data on which the base percentiles are calculated.

  The \code{n} argument specifies the size of the window used when
  computing the percentiles used in \code{\link{climdex.tx10p}},
  \code{\link{climdex.tn10p}}, \code{\link{climdex.tx90p}}, and
  \code{\link{climdex.tn90p}}.

  The \code{northern.hemisphere} argument specifies whether the data
  came from the northern hemisphere. If FALSE, data is assumed to have
  come from the southern hemisphere. This is used when computing growing
  season length; if the data is from the southern hemisphere, growing
  season length is the growing season starting in the beginning of July
  of the year indicated, running to the end of June of the following
  year.

  The \code{pad.data.with.first.last.values} argument specifies whether
  to pad the data passed into the baseline quantile routine with the
  first and last values. If TRUE, the first (at the beginning of the
  series) and last (at the end of the series) values will be used to pad
  the beginning and ends of this series. If FALSE, either NA or the
  values for the previous two (at the beginning) and last two (at the
  end) days of data will be used.
  
  The \code{temp.quantiles.notbase} argument allows the user to supply
  pre-computed quantiles for data outside the base period. This is a
  list consisting of four named vectors of length 365 (or 360, in the
  case of 360 day calendars). These vectors contain the thresholds for
  each day for the 10th percentile of tmax and tmin (tx10thresh and
  tn10thresh) and for the 90th percentile of tmax and tmin (tx90thresh
  and tn90thresh).
  
  The \code{prec.quantiles} argument allows the user to specify
  pre-computed precipitation quantiles. This is a vector containing two
  named elements -- r95thresh (the 95th percentile of wet days) and
  r99thresh (the 99th percentile of wet days).

  The \code{temp.quantiles.base} argument allows the user to specify
  pre-computed quantiles for data inside the base period. This is a list
  consisting of four named vectors of dimensions 365 (or 360, in the
  case of 360 day calendars) by the number of base years by the number
  of base years less one; two for the 10th percentile of tmax and tmin
  (tx10thresh and tn10thresh) and two for the 90th percentile of tmax
  and tmin (tx90thresh and tn90thresh). Each of these vectors contain a
  set of thresholds for each day for each year; this set contains
  quantiles computed with the current year replaced with one of the
  other years of data.
}

\examples{
## Create a climdexInput object from some data already loaded in and
## ready to go.

## Parse the dates into PCICt.
tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
"jday")]), format="\%Y \%j", cal="gregorian")
tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
"jday")]), format="\%Y \%j", cal="gregorian")
prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
"jday")]), format="\%Y \%j", cal="gregorian")

## Load the data in.
ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))

## Alternatively, one could load the data in directly using
## climdexInput.csv from CSV files. This would create a climdexInput
## object from a set of filenames (already stored as variables), with a
## different date format.
\donttest{ci.csv <- climdexInput.csv(tmax.filename, tmin.filename,
prec.filename, date.types=list(list(fields=c("date"), format="\%Y-\%m-\%d")))}
}

\seealso{
  \code{\link{climdex.pcic-package}}, \code{\link{strptime}}.
}

\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
}  

\keyword{ts}
\keyword{climate}