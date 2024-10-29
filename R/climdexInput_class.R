#' @importFrom utils head read.csv
#' 

## Class definition declaration
#' @title climdexInput
#' 
#' @description
#' The climdexInput class contains all the data necessary to compute the
#' climdex indices.
#' 
#' @details
#' The \code{climdexInput} class consists of all the data necessary to compute
#' the climdex indices. Users will not need to modify any of the slots in this
#' class. That being said, users may want or need to repurpose this data for
#' further analysis. The following description of the data is aimed at that
#' audience.
#' 
#' The \code{data} slot contains time series' of daily data of equal length for
#' each of the provided variables. Missing days have been replaced with NA.
#' The \code{dates} slot is the corresponding series of dates (of type PCICt)
#' for the daily data.
#' 
#' The \code{quantiles} slot contains quantiles used for computing the
#' tn/tx 10/90p indices, w/csdi, r95ptot, and r99ptot. If precipitation data
#' is supplied, the 'prec' member contains the 95th and 99th percentile values
#' for precipitation within the base period. For tmin and tmax, if present each
#' will have a corresponding member in the slot. Within each of these, there
#' will be an 'inbase' and 'outbase' member, corresponding to thresholds to be
#' used within the base period (inbase) and outside the base period (outbase).
#' The 'inbase' member consists of one percentile for each day of the year,
#' computed using an n-day (default is 5-day) running window surrounding that
#' day. These percentiles are computed for at least the 10th and 90th
#' percentile of the data. For the 'outbase' member, given n years
#' of data to use as the base period, there are n * (n - 1) sets of daily
#' quantiles of the same type as those in 'inbase'.
#' 
#' To ease computation of monthly, seasonal and annual data, \code{date.factors} 
#' contains date factors which group data into annual, seasonal and monthly time
#' buckets. They are of the same length as the time series and can be reused
#' for computation of any annual, seasonal or monthly aggregates.
#' 
#' The climdexInput class also includes NA masks for seasonal, monthly
#' and annual as parts of the \code{namasks} slot. Each of these masks consist
#' of a vector of numbers of the same length as the monthly, seasonal or annual output
#' data. The values used are 1 to signify that the data meets the QC criteria,
#' and NA to signify it does not. Years with more than (by default) 15 days
#' missing, months with more than (by default) 3 days missing, 
#' and seasons with more than (by default) 6 days missing are
#' considered to be of poor quality and are masked here with NA. These
#' thresholds can be set when instantiating the object, and are stored in the
#' \code{max.missing.days} slot.
#' 
#' Seasons are considered valid only when all three constituent months are valid,
#' and the sum of NA days for the season is no greater than the max.missing.days
#' seasonal argument. This means that a time series starting in January and ending
#' in December will have NA winter seasons at both ends of the series.
#'
#' Seasons are defined as the meteorological seasons:
#'
#' - 'winter': December, January, February
#' - 'spring': March, April, May
#' - 'summer': June, July, August
#' - 'autumn': September, October, November
#' 
#' The \code{base.range} slot contains vector of type PCICt containing the
#' first and last day included in the baseline.
#' 
#' The \code{northern.hemisphere} slot contains a boolean indicating whether
#' the data came from the northern hemisphere. If FALSE, data is assumed to
#' have come from the southern hemisphere. This is used when computing growing
#' season length; if the data is from the southern hemisphere, growing season
#' length is the growing season starting in the beginning of July of the year
#' indicated, running to the end of June of the following year.
#' 
#' The \code{max.missing.days} slot is a vector consisting of 'annual'
#' (the number of days that can be missing in a year), 'monthly' (the
#' number of days that can be missing in a month and 'seasonal' (the
#' number of days that can be missing in a season. If one month in a year fails
#' the test, the corresponding year will be omitted.
#' 
#' @name climdexInput
#' @aliases climdexInput-class
#' @docType class
#' @section Slots: \describe{
#' \item{data}{Time series of supplied data variables.}
#' \item{quantiles}{Threshold quantiles used for threshold-based indices.}
#' \item{namasks}{Data quality masks for annual, seasonal and monthly data.}
#' \item{dates}{Date sequence (type PCICt) corresponding to temperature and
#' precipitation data.}
#' \item{jdays}{Julian days for the date sequence.}
#' \item{base.range}{Date range (type PCICt) of baseline period.}
#' \item{date.factors}{Factors used for creation of annual, seasonal and monthly indices.}
#' \item{northern.hemisphere}{Boolean used when computing growing season
#' length.}
#' \item{max.missing.days}{Maximum number of missing days of data for annual, seasonal
#' and monthly data.}
#' }
#' @seealso \code{\link{climdexInput.csv}}, \code{\link{climdexInput.raw}}.
#' @keywords climate ts
#' @examples
#' library(PCICt)
#' 
#' ## Parse the dates into PCICt.
#' tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' 
#' ## Load the data in.
#' ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#' 
#' @export
setClass("climdexInput",
         representation(data = "list",
                        quantiles = "environment",
                        namasks = "list",
                        dates = "PCICt",
                        jdays = "numeric",
                        base.range = "PCICt",
                        date.factors = "list",
                        northern.hemisphere = "logical",
                        max.missing.days = "numeric"),
         validity=valid.climdexInput
)

#' @title Method for creating climdexInput object from vectors of data
#' 
#' @description
#' This function creates a climdexInput object from data already ingested into
#' R.
#' 
#' @details
#' This function takes input climate data at daily resolution, and produces as
#' output a ClimdexInput data structure. This data structure can then be passed
#' to any of the routines used to compute the Climdex indices. The indices
#' themselves are specified on the webpage cited in the references section.
#' The \code{base.range} argument is a pair of 4 digit years which bound the
#' data on which the base percentiles are calculated.
#' 
#' The \code{tmax}, \code{tmin}, and \code{prec} arguments are numeric vectors
#' containing the data on which the indices are to be computed. The units are
#' assumed to be degrees C for temperature, and mm/day for precipitation.
#'
#' The \code{tmax.dates}, \code{tmin.dates}, and \code{prec.dates} arguments
#' are vectors of type \code{PCICt}.
#' 
#' The \code{n} argument specifies the size of the window used when computing
#' the percentiles used in \code{\link{climdex.tx10p}},
#' \code{\link{climdex.tn10p}}, \code{\link{climdex.tx90p}}, and
#' \code{\link{climdex.tn90p}}.
#' 
#' The \code{northern.hemisphere} argument specifies whether the data came from
#' the northern hemisphere. If FALSE, data is assumed to have come from the
#' southern hemisphere. This is used when computing growing season length; if
#' the data is from the southern hemisphere, growing season length is the
#' growing season starting in the beginning of July of the year indicated,
#' running to the end of June of the following year.
#' 
#' The \code{quantiles} argument allows the user to supply pre-computed quantiles.
#' This is a list consisting of quantiles for each variable.
#' 
#' For each temperature variable, there are separate lists of quantiles for 
#' inbase and outbase, with these names. In both cases, quantiles within these
#' lists are named q10 for the 10th percentile and q90 for the 90th percentile.
#' Other percentiles would be named qnn for the nnth percentile. For the
#' outbase quantiles, each element in the list is a vector of length 365 (or 360
#' in the case of 360-day calendars), corresponding to one value for each day of
#' the year. For the inbase quantiles, each element in the list is an array of
#' dimensions \code{[365 or 360, nyr, nyr - 1]}, where nyr is the number of years in
#' the base period. Each value corresponds to a quantile for each day, for each
#' year, with a particular year replaced.
#'
#' For precipitation variables, there is a named vector of quantiles, consisting
#' of at least q95 and q99. 
#'
#' The \code{temp.qtiles} and \code{prec.qtiles} arguments allow the user to
#' modify the quantiles calculated. For example, specifying
#' temp.qtiles=c(0.10, 0.50, 0.90) would calculate the 10th, 50th, and 90th
#' percentiles for temperature.
#'
#' The \code{min.base.fraction.present} argument specifies the minimum fraction
#' of data which must be present for a quantile to be calculated for a 
#' particular day. If the fraction of data present is less than this threshold, 
#' the quantile for that day will be set to NA.
#'
#' The \code{max.missing.days} argument is a vector consisting of 'annual'
#' (the number of days that can be missing in a year), 'monthly' (the
#' number of days that can be missing in a month and 'seasonal' (the
#' number of days that can be missing in a season. If one month in a year fails
#' the test, the corresponding year will be omitted.
#' 
#' @seealso \code{\link{climdex.pcic-package}}, \code{\link{strptime}}.
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @keywords ts climate
#'
#' @template climdexInput_raw_help1 
#' @template climdexInput_raw_params
#' @template climdexInput_common_params
#' @param northern.hemisphere Whether this point is in the northern hemisphere.
#' @param quantiles Threshold quantiles for supplied variables.
#' @param max.missing.days Vector containing thresholds for number of days
#' allowed missing per year (annual), per month (monthly) and per season (seasonal).
#' @return An object of class \code{\link{climdexInput-class}} for use with
#' other climdex methods.
#' @note Units are assumed to be mm/day for precipitation and degrees Celsius
#' for temperature. No units conversion is performed internally.
#' 
#' @examples
#' library(PCICt)
#'
#' ## Create a climdexInput object from some data already loaded in and
#' ## ready to go.
#' 
#' ## Parse the dates into PCICt.
#' tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' 
#' ## Load the data in.
#' ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#'
#' @export 
climdexInput.raw <- function(tmax=NULL, tmin=NULL, prec=NULL, tmax.dates=NULL, tmin.dates=NULL, prec.dates=NULL,
                             base.range=c(1961, 1990), n=5, northern.hemisphere=TRUE,
                             tavg=NULL, tavg.dates=NULL, quantiles=NULL, temp.qtiles=c(0.10, 0.90), prec.qtiles=c(0.95, 0.99), max.missing.days=c(annual=15, monthly=3, seasonal=6), min.base.data.fraction.present=0.1) {
  ## Make sure all of these arguments are valid...
  check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n, tavg, tavg.dates)
  
  stopifnot(length(max.missing.days) == 3 && all(c("annual", "monthly", "seasonal") %in% names(max.missing.days)))
  stopifnot(is.numeric(min.base.data.fraction.present) && length(min.base.data.fraction.present) == 1)
  
  d.list <- list(tmin.dates, tmax.dates, prec.dates, tavg.dates)
  all.dates <- do.call(c, d.list[!sapply(d.list, is.null)])
  last.day.of.year <- get.last.monthday.of.year(all.dates)
  cal <- attr(all.dates, "cal")
  
  ## Convert base range (in years) to PCICt
  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)
  bs.date.series <- seq(bs.date.range[1], bs.date.range[2], by="day")
  
  ## Get dates for normal data
  new.date.range <- as.PCICt(paste(as.numeric(format(range(all.dates), "%Y", tz="GMT")), c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")
  jdays <- get.jdays.replaced.feb29(get.jdays(date.series))
  season_with_year <- classify_meteorological_season_with_year(date.series)
  
  ## Factors for dividing data up
  date.factors <- list(
    annual = factor(format(date.series, format = "%Y", tz = "GMT")),
    monthly = factor(format(date.series, format = "%Y-%m", tz = "GMT")),
    seasonal = factor(season_with_year, levels = unique(season_with_year))
  )
  ## Filled data...
  var.list <- c("tmax", "tmin", "prec", "tavg")
  present.var.list <- var.list[sapply(var.list, function(x) !is.null(get(x)))]
  filled.list <- sapply(present.var.list, function(x) { return(create.filled.series(get(x), trunc(get(paste(x, "dates", sep="."))), date.series)) }, simplify=FALSE)
  if(is.null(tavg) && !is.null(tmin) && !is.null(tmax))
    filled.list$tavg <- (filled.list$tmax + filled.list$tmin) / 2
  
  ## Establish some truth values for later use in logic...
  days.threshold <- 359
  present.dates <- sapply(present.var.list, function(x) get(paste(x, "dates", sep=".")))
  quantile.dates <- list(tmax=tmax.dates, tmin=tmin.dates, prec=prec.dates)
  days.in.base <- sapply(quantile.dates, get.num.days.in.range, bs.date.range)
  
  ## Check that provided quantiles, if any, are valid
  check.quantile.validity(quantiles, present.var.list, days.in.base)
  
  data.in.base.period <- any(days.in.base != 0)
  have.quantiles <- all(present.var.list %in% names(quantiles))
  
  ## NA masks
  namasks <- list(
    annual = lapply(filled.list, get.na.mask, date.factors$annual, max.missing.days["annual"]),
    monthly = lapply(filled.list, get.na.mask, date.factors$monthly, max.missing.days["monthly"]),
    seasonal = lapply(filled.list, get.na.mask, date.factors$seasonal, max.missing.days["seasonal"]))
  namasks$annual <- lapply(names(namasks$annual), function(v) {
    d <- namasks$annual[[v]] * as.numeric(tapply(namasks$monthly[[v]], rep(seq_along(namasks$annual[[v]]), each = 12), prod))
    dimnames(d) <- dim(d) <- NULL
    d
  })
  names(namasks$annual) <- names(namasks$seasonal) <- names(namasks$monthly)
  
  
  season_month_counts <- sapply(unique(date.factors$seasonal), function(season) {
    length(unique(date.factors$monthly[date.factors$seasonal == season]))
  })
  data.vars <- present.var.list
  if (!is.null(namasks$seasonal$tavg)){
    data.vars <- c(data.vars, "tavg")
  }
  for (var in data.vars) {
    seasonal_namasks <- namasks$seasonal[[var]]
    na_months <- unique(date.factors$monthly)[is.na(namasks$monthly[[var]])]
    seasons_of_na_months <- unique(date.factors$seasonal[date.factors$monthly %in% na_months])
    seasonal_namasks[unique(date.factors$seasonal) %in% seasons_of_na_months] <- NA
    # Identify and set NA for seasons with less than 3 months
    for (season in seq_along(season_month_counts) ) {
      if (!is.na(season_month_counts[season]) && season_month_counts[season] < 3) {
        seasonal_namasks[season] <- NA
      }
    }
    namasks$seasonal[[var]] <- seasonal_namasks
  } 
  
  
  ## Pad data passed as base if we're missing endpoints...
  if(!have.quantiles) {
    quantiles <- new.env(parent=emptyenv())
    
    if(days.in.base['tmax'] > days.threshold)
      delayedAssign("tmax", get.temp.var.quantiles(filled.list$tmax, date.series, bs.date.series, temp.qtiles, bs.date.range, n, TRUE, min.base.data.fraction.present), assign.env=quantiles)
    if(days.in.base['tmin'] > days.threshold)
      delayedAssign("tmin", get.temp.var.quantiles(filled.list$tmin, date.series, bs.date.series, temp.qtiles, bs.date.range, n, TRUE, min.base.data.fraction.present), assign.env=quantiles)
    if(days.in.base['prec'] > days.threshold)
      delayedAssign("prec", get.prec.var.quantiles(filled.list$prec, date.series, bs.date.range, prec.qtiles), assign.env=quantiles)
  } else {
    quantiles <- as.environment(quantiles)
  }
  
  return(new("climdexInput", data=filled.list, quantiles=quantiles, namasks=namasks, dates=date.series, jdays=jdays, base.range=bs.date.range, date.factors=date.factors, northern.hemisphere=northern.hemisphere, max.missing.days=max.missing.days))
}

#' @title Method for creating climdexInput object from CSV files
#' 
#' @description
#' This function creates a climdexInput object from data in CSV files.
#' 
#' @details
#' This function takes input climate data in CSV files at daily resolution,
#' and produces as output a ClimdexInput data structure. This data structure
#' can then be passed to any of the routines used to compute the Climdex
#' indices. The indices themselves are specified on the webpage cited in the
#' references section.
#'
#' Any of tmin.file (daily minimum temperature), tmax.file (daily maximum
#' temperature), tavg.file (daily mean temperature), and prec.file (daily
#' precipitation) can be passed in. tavg will be derived from the mean of
#' tmax and tmin if it is not supplied. If any of tmin.file, tmax.file, and
#' prec.file are not supplied, the set of indices which can be calculated will
#' be limited to indices which do not involve the missing variables.
#' 
#' The \code{tmax.file}, \code{tmin.file}, and \code{prec.file} arguments
#' should be names of CSV files containing dates and the data on which the
#' indices are to be computed. The units are assumed to be degrees C for
#' temperature, and mm/day for precipitation.
#' 
#' The \code{data.columns} argument is a vector consisting of named items tmax,
#' tmin, and prec. These named items are used as the column names in their
#' respective files when loading in CSV.
#' 
#' The \code{cal} argument is a textual description of the calendar type, as
#' described in the documentation for \code{\link{as.PCICt}}.
#' 
#' The \code{date.types} argument is a list of lists containing two named
#' items: \code{fields}, and \code{format}. The \code{fields} item is a vector
#' of names consisting of the columns to be concatenated together with spaces.
#' The \code{format} item is a date format as taken by \code{strptime}.
#'
#' For more details on arguments, see \code{\link{climdexInput.raw}}.
#'
#' @seealso \code{\link{climdex.pcic-package}}, \code{\link{climdexInput.raw}}.
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @keywords ts climate
#' 
#' @param tmax.file Name of file containing daily maximum temperature data.
#' @param tmin.file Name of file containing daily minimum temperature data.
#' @param prec.file Name of file containing daily total precipitation data.
#' @param tavg.file Name of file containing daily mean temperature data.
#' @param data.columns Column names for tmin, tmax, and prec data.
#' @param date.types Column names for tmin, tmax, and prec data (see notes).
#' @param na.strings Strings used for NA values; passed to
#' \code{\link{read.csv}}.
#' @param cal The calendar type used in the input files.
#' @template climdexInput_common_params
#' @param northern.hemisphere Whether this point is in the northern hemisphere.
#' @param quantiles Threshold quantiles for supplied variables.
#' @param max.missing.days Vector containing thresholds for number of days
#' allowed missing per year (annual), per month (monthly) and per season (seasonal).
#' @return An object of class \code{\link{climdexInput-class}} for use with
#' other climdex methods.
#' @note Units are assumed to be mm/day for precipitation and degrees Celsius
#' for temperature. No units conversion is performed internally.
#' @examples
#' ## This would create a climdexInput object from a set of filenames (already
#' ## stored as variables), with a different date format.
#' \dontrun{ci.csv <- climdexInput.csv(tmax.filename, tmin.filename,
#' prec.filename, date.types=list(list(fields=c("date"), format="%Y-%m-%d")))}
#'
#' @export
climdexInput.csv <- function(tmax.file=NULL, tmin.file=NULL, prec.file=NULL,
                             data.columns=list(tmin="tmin", tmax="tmax", prec="prec"), base.range=c(1961, 1990),
                             na.strings=NULL, cal="gregorian", date.types=NULL, n=5, northern.hemisphere=TRUE,
                             tavg.file=NULL, quantiles=NULL, temp.qtiles=c(0.10, 0.90), prec.qtiles=c(0.95, 0.99), max.missing.days=c(annual=15, monthly=3, seasonal=6), min.base.data.fraction.present=0.1) {
  get.and.check.data <- function(fn, datacol) {
    if(!is.null(fn)) {
      dat <- read.csv(fn, na.strings=na.strings)
      if(!(datacol %in% names(dat)))
        stop("Data column not found in tmin data.")
      return(list(dat=dat[!is.na(dat[,datacol]),datacol],  dates=get.date.field(dat, cal, date.types)))
    }
    return(list(dat=NULL, dates=NULL))
  }
  if(missing(date.types))
    date.types <- list(list(fields=c("year", "jday"), format="%Y %j"),
                       list(fields=c("year", "month", "day"), format="%Y %m %d"))
  else
    if(any(!sapply(date.types, function(x) { return(sum(c("fields", "format") %in% names(x)) == 2 && is.character(x$fields) && is.character(x$format)) } )))
      stop("Invalid date.types specified. See ?climdexInput.csv .")
  
  tmin <- get.and.check.data(tmin.file, data.columns$tmin)
  tmax <- get.and.check.data(tmax.file, data.columns$tmax)
  tavg <- get.and.check.data(tavg.file, data.columns$tavg)
  prec <- get.and.check.data(prec.file, data.columns$prec)
  
  return(climdexInput.raw(tmax=tmax$dat, tmin=tmin$dat, prec=prec$dat, tmax.dates=tmax$dates, tmin.dates=tmin$dates, prec.dates=prec$dates, base.range=base.range, n=n, northern.hemisphere=northern.hemisphere, tavg=tavg$dat, tavg.dates=tavg$dates, quantiles=quantiles, temp.qtiles=temp.qtiles, prec.qtiles=prec.qtiles, max.missing.days=max.missing.days, min.base.data.fraction.present=min.base.data.fraction.present))
}
