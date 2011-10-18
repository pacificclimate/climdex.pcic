\name{climdexInput.raw, climdexInput.csv}
\alias{climdexInput.csv}
\alias{climdexInput.raw}
\title{Methods for Creating climdexInput Objects}
\description{
  These functions create climdexInput objects for use with the climdex methods.
}
\usage{
climdexInput.raw(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates,
  base.range=c(1961, 1990), n=5)

climdexInput.csv(tmax.file, tmin.file, prec.file,
  data.columns=list(tmin="tmin", tmax="tmax", prec="prec"),
  base.range=c(1961, 1990), na.strings=NULL, cal="gregorian",
  date.types=NULL, n=5)
}
\arguments{
  \item{tmax}{Daily maximum temperature data.}
  \item{tmin}{Daily minimum temperature data.}
  \item{prec}{Daily total precipitation data.}
  \item{tmax.dates}{Dates for the daily maximum temperature data.}
  \item{tmin.dates}{Dates for the daily minimum temperature data.}
  \item{prec.dates}{Dates for the daily total precipitation data.}
  \item{tmax.file}{Name of file containing daily maximum temperature data.}
  \item{tmin.file}{Name of file containing daily minimum temperature data.}
  \item{prec.file}{Name of file containing daily total precipitation data.}
  \item{data.columns}{Column names for tmin, tmax, and prec data.}
  \item{date.types}{Column names for tmin, tmax, and prec data (see notes).}
  \item{na.strings}{Strings used for NA values; passed to \code{\link{read.csv}}.}
  \item{cal}{The calendar type used in the input files.}
  \item{base.range}{Years to use for the baseline.}
  \item{n}{Number of days to use as window for daily quantiles.}
}

\details{
  These functions take input climate data, and produce as output a data
  structure which can be used with the other methods in this package to
  compute all of the Climdex indices as specified on the webpage cited
  in the references section.

  \code{climdexInput.csv} should be used when your data is already in a
  CSV file. \code{climdexInput.raw} should be used in all other circumstances, as
  it is much more flexible.
  
  \code{climdexInput.csv} takes as arguments three input filenames, one for
  each of daily maximum temperature, minimum temperature, and total
  preciptation. Optionally, it also takes a list containing the names of
  the data columns for each of the variables stored in those files, a
  vector of NA strings for missing value handling, the baseline range,
  the window size, the calendar type, and a data structure containing
  the information needed to parse the dates out of the file.

  \code{climdexInput.raw} takes as arguments the input daily maximum
  temperature, minimum temperature, and total precipitation data along
  with the associated series of dates and optionally the baseline range
  and window size.
}

\value{
  An object of class \code{\link{climdexInput-class}} for use with other climdex methods.
}

\note{
  Units are assumed to be mm/day for precipitation and degrees Celsius
  for temperature. No units conversion happens internally.

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