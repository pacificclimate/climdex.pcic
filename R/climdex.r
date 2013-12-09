library(caTools)
library(PCICt)

#' Get the last month and day of the year
#'
#' Get the last month and day of the year as a character sting, separated by
#' the specified separator.
#'
#' This is a utility function necessitated by 360-day calendars. Works on PCICt objects.
#'
#' @param d An exemplar date.
#' @param sep Separator to use.
#' @return A string (like "12-30", or "12-31")
#' 
#' @examples
#' last.mday <- get.last.monthday.of.year(as.PCICt("2011-01-01", cal="360"))
#' 
#' @export
get.last.monthday.of.year <- function(d, sep="-") {
  if(!is.null(attr(d, "months"))) paste("12", attr(d, "months")[12], sep=sep) else paste("12", "31", sep=sep)
}

## Lower overhead version of tapply
tapply.fast <- function (X, INDEX, FUN = NULL, ..., simplify = TRUE) {
  FUN <- if (!is.null(FUN))
    match.fun(FUN)
  
  if(!is.factor(INDEX))
    stop("INDEX must be a factor.")
  
  if (length(INDEX) != length(X))
    stop("arguments must have same length")
  
  if (is.null(FUN))
    return(INDEX)
  
  namelist <- levels(INDEX)
  ans <- lapply(split(X, INDEX), FUN, ...)
  
  ans <- unlist(ans, recursive = FALSE)
  names(ans) <- levels(INDEX)
  return(ans)
}

## Check that climdexInput data structure is valid.
valid.climdexInput <- function(x) {
  temp.quantiles <- c(10, 90)
  prec.quantiles <- c(95, 99)
  errors <- c()

  separate.base <- c(tmax=T, tmin=T, tavg=T, prec=F)
  present.data.vars <- names(x@data)
  length.check.slots <- c("dates", "jdays", "annual.factor", "monthly.factor")
  data.lengths <- c(sapply(x@data, length), sapply(length.check.slots, function(y) length(slot(x, y))))
  quantiles <- list(tmax=temp.quantiles, tmin=temp.quantiles, prec=prec.quantiles)
  
  if(!all(data.lengths == max(data.lengths)))
    errors <- c(errors, "Data fields, dates, and factors must all be of the same length")

  ## Check that namask.mon and namask.ann have columns for each of the variables
  if(!all(present.data.vars %in% names(x@namask.ann) & present.data.vars %in% names(x@namask.mon)))
    errors <- c(errors, "NA mask for monthly and annual must contain data for all variables supplied.")

  ## Check that appropriate thresholds are present.
  need.base.data <- get.num.days.in.range(x@dates, x@base.range) > 0
  errors <- do.call(c, c(list(errors), lapply(intersect(present.data.vars, c("tmax", "tmin", "prec")), function(n) {
    if(is.null(quantiles[n]))
      return(NULL)
    if(!(n %in% ls(envir=x@quantiles)))
      return(paste("Quantiles for", n, "are missing.", sep=""))
    return(NULL)
  })))

  if(length(x@northern.hemisphere) != 1)
    errors <- c(errors, "northern.hemisphere must be of length 1.")
  
  if(length(errors) == 0)
    return(TRUE)
  else
    return(errors)
}

## Class definition declaration
#' climdexInput
#' 
#' The climdexInput class contains all the data necessary to compute the
#' climdex indices.
#' 
#' The \code{climdexInput} class consists of all the data necessary to compute
#' the climdex indices. Users will not need to modify any of the slots in this
#' class. That being said, users may want or need to repurpose this data for
#' further analysis. The following description of the data is aimed at that
#' audience.
#' 
#' The \code{tmax}, \code{tmin}, \code{tavg}, and \code{prec} slots are time
#' series of daily data of equal length and without any missing days, with NAs
#' in place of data where no data was present. The \code{dates} slot is the
#' corresponding series of dates (of type PCICt) for the daily data.
#' 
#' To ease computation of monthly and annual data, \code{monthly.factor} and
#' \code{annual.factor} are slots in the data structure. They are also of the
#' same length as the time series. These can be reused for computation of any
#' annual or monthly aggregates.
#' 
#' The climdexInput class also includes NA masks for both monthly
#' (\code{namask.mon}) and annual (\code{namask.ann}) data. These masks consist
#' of a vector of numbers of the same length as the monthly or annual output
#' data. The values used are 1 to signify that the data meets the QC criteria,
#' and NA to signify it does not. Years with more than 15 days missing, and
#' months with more than 3 days missing, are considered to be of poor quality
#' and are masked here with NA.
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
#' The \code{pctile} slot contains the 95th and 99th percentile values for the
#' precipitation percentiles, computed on the base period.
#' 
#' The \code{running.pctile.notbase} slot contains the data necessary for
#' computing temperature based percentiles outside of the base period. This
#' consists of one percentile for each day of the year, computed using an n-day
#' (default is 5-day) running window surrounding that day. These percentiles
#' are computed for both the 10th and 90th percentile for daily maximum and
#' minimum temperature.
#' 
#' The \code{running.pctile.base} slot contains the data necessary for
#' computing temperature based percentiles inside the base period. As this is a
#' somewhat unpleasant operation, so are the data requirements. Given n years
#' of data to use as the base period, there are n * (n - 1) sets of daily
#' quantiles of the same type as those for \code{running.pctile.notbase}.
#' 
#' @name climdexInput
#' @aliases climdexInput-class
#' @docType class
#' @section Slots: \describe{ \item{data}{Time series of supplied data
#' variables.} \item{quantiles}{Threshold quantiles used for
#' threshold-based indices.}\item{namask.ann}{Data quality mask for
#' annual data.} \item{namask.mon}{Data quality mask for monthly data.}
#' \item{dates}{Date sequence (type PCICt) corresponding to temperature and
#' precipitation data.}
#' \item{jdays}{Julian days for the date sequence.}
#' \item{base.range}{Date range (type PCICt) of baseline period.}
#' \item{annual.factor}{Factor used for creation of annual indices.}
#' \item{monthly.factor}{Factor used for creation of monthly indices.}
#' \item{northern.hemisphere}{Boolean used when computing growing season
#' length.} }
#' @seealso \code{\link{climdexInput.csv}}, \code{\link{climdexInput.raw}}.
#' @keywords climate ts
#' @examples
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
                        namask.ann = "list",
                        namask.mon = "list",
                        dates = "PCICt",
                        jdays = "numeric",
                        base.range = "PCICt",
                        annual.factor = "factor",
                        monthly.factor = "factor",
                        northern.hemisphere = "logical"),
         validity=valid.climdexInput
         )

## Returns PCICt field or dies
get.date.field <- function(input.data, cal, date.types) {
  valid.date.types <- sapply(date.types, function(x) { return(!inherits(try(input.data[,x$fields], silent=TRUE), "try-error")) })

  if(sum(valid.date.types) == 0) {
    stop("Could not find a workable set of date fields")
  }

  date.type <- date.types[[which(valid.date.types)[1]]]
  date.strings <- do.call(paste, input.data[,date.type$fields])
  return(as.PCICt(date.strings, format=date.type$format, cal=cal))
}

## Creates a filled series given the data, dates, and new date sequence to be used.
create.filled.series <- function(data, data.dates, new.date.sequence) {
  new.data <- rep(NA, length(new.date.sequence))
  data.in.new.data <- (data.dates >= new.date.sequence[1]) & (data.dates <= new.date.sequence[length(new.date.sequence)])
  indices <- floor(as.numeric(data.dates[data.in.new.data] - new.date.sequence[1], units="days")) + 1
  new.data[indices] <- data[data.in.new.data]
  return(new.data)
}

## Get julian day of year
get.jdays <- function(dates) {
  return(as.POSIXlt(dates)$yday + 1)
}

## Get year
get.years <- function(dates) {
  return(as.POSIXlt(dates)$year + 1900)
}

## Get month number
get.months <- function(dates) {
  return(as.POSIXlt(dates)$mon + 1)
}

## Juggle the list so that day 366 == day 365
get.jdays.replaced.feb29 <- function(jdays) {
  indices <- which(jdays == 366)
  if(length(indices) > 0)
    jdays[rep(indices, each=366) + -365:0] <- c(1:59, 59, 60:365)
  jdays
}

## Get set of days for bootstrap use
get.bootstrap.set <- function(dates, bootstrap.range, win.size) {
  bootstrap.win.range <- get.bootstrap.windowed.range(bootstrap.range, win.size)
  dpy <- ifelse(is.null(attr(dates, "dpy")), 365, attr(dates, "dpy"))
  return(dates >= bootstrap.win.range[1] & dates <= bootstrap.win.range[2] & (dpy == 360 | format(dates, format="%m-%d", tz="GMT") != "02-29"))
}

## Get bootstrap date range
## Input is: bootstrap range, vector of 2 PCICt, win.size is window size in days (single integer)
get.bootstrap.windowed.range <- function(bootstrap.range, win.size) {
  ## Changed due to a bug in PCICt
  ##window <- as.difftime(floor(win.size / 2), units="days")
  window <- floor(win.size / 2) * 86400
  return(c(bootstrap.range[1] - window, bootstrap.range[2] + window))
}

## Calculate a running quantile on the data set over the bootstrap range.
## If get.bootstrap.data is TRUE, use the Zhang boostrapping method described in Xuebin Zhang et al's 2005 paper, "Avoiding Inhomogeneity in Percentile-Based Indices of Temperature Extremes" J.Clim vol 18 pp.1647-1648, "Removing the 'jump'".
## Expects PCICt for all dates
zhang.running.qtile <- function(x, dates.base, qtiles, bootstrap.range, include.mask=NULL, n=5, pad.data.with.first.last.values=FALSE, get.bootstrap.data=FALSE) {
  inset <- get.bootstrap.set(dates.base, bootstrap.range, n)
  dpy <- ifelse(is.null(attr(dates.base, "dpy")), 365, attr(dates.base, "dpy"))
  nyears <- floor(sum(inset) / dpy)
  
  if(!is.null(include.mask))
    x[include.mask] <- NA

  bs.data <- x[inset]
  window <- floor(n / 2)
  if(pad.data.with.first.last.values) {
    bs.data[1:window] <- bs.data[window + 1]
    bs.data[length(bs.data) - 0:(window - 1)] <- bs.data[length(bs.data) - window]
  }

  qdat <- NULL
  if(get.bootstrap.data) {
    d <- .Call("running_quantile_windowed_bootstrap", bs.data, n, qtiles, dpy, DUP=FALSE, PACKAGE='climdex.pcic')
    dim(d) <- c(dpy, nyears, nyears - 1, length(qtiles))
    qdat <- lapply(1:length(qtiles), function(x) { d[,,,x] })
  } else {
    res <- running.quantile(bs.data, n, qtiles, dpy)
    qdat <- lapply(1:length(qtiles), function(x) { res[,x] })
  }
  names(qdat) <- paste("q", qtiles * 100, sep="")
  return(qdat)
}

## Get NA mask given threshold and split factor
get.na.mask <- function(x, f, threshold) {
  return(c(1, NA)[1 + as.numeric(tapply.fast(is.na(x), f, function(y) { return(sum(y) > threshold) } ))])
}

## Get number of days within range
get.num.days.in.range <- function(x, date.range) {
  return(sum(x >= date.range[1] & x <= date.range[2]))  
}


## Check that arguments to climdexInput.raw et al are complete enough and valid enough.
check.basic.argument.validity <- function(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range=c(1961, 1990), n=5, tavg=NULL, tavg.dates=NULL) {
  check.var <- function(var, var.dates, var.name) {
    if(is.null(var) != is.null(var.dates))
      stop(paste("If passing in", var, ", must pass in", var, "dates too.."))
    if(!is.null(var.dates) && length(var) != length(var.dates))
      stop(paste("Length of", var.name, "data and dates do not match."))
    if(!is.null(var.dates) && !inherits(var.dates, "PCICt"))
      stop(paste(var.name, "dates must be of class PCICt."))
    if(!is.null(var) && !is.numeric(var))
      stop(paste(var.name, "must be of type numeric."))
  }
    
  check.var(tmax, tmax.dates, "tmax")
  check.var(tmin, tmin.dates, "tmin")
  check.var(tavg, tavg.dates, "tavg")
  check.var(prec, prec.dates, "prec")

  if(sum(c(is.null(tmax), is.null(tmin), is.null(prec), is.null(tavg))) == 0)
    stop("Must supply at least one variable to calculate indices upon.")

  if(!(length(base.range) == 2 && is.numeric(base.range)))
    stop("Invalid base date range; expecting vector of 2 numeric years.")

  if(!is.numeric(n) || length(n) != 1)
    stop("n must be numeric and of length 1.")
  
  if(n != 5)
    warning("Use of n != 5 varies from the Climdex definition. Use at your own risk.")
}

## Check validity of quantile input.
check.quantile.validity <- function(quantiles, present.vars, days.in.base) {
  if(is.null(quantiles))
    return()
  
  if(class(quantiles) != "list")
    stop("Provided quantiles must be a list.")
  
  if(!all(present.vars %in% names(quantiles)))
    stop("Quantiles must be present for all variables provided.\n")

  if(!all(sapply(quantiles[names(quantiles) %in% intersect(present.vars, c("tmax", "tmin"))], function(x) { "outbase" %in% names(x) && all(c("q10", "q90") %in% names(x$outbase)) })))
    stop("Temperature out-of-base quantiles must contain 10th and 90th percentiles.\n")

  if(any(days.in.base > 0) && !all(sapply(quantiles[names(quantiles) %in% intersect(intersect(present.vars, c("tmax", "tmin")), names(days.in.base)[days.in.base > 0])], function(x) { "inbase" %in% names(x) && all(c("q10", "q90") %in% names(x$inbase)) })))
    stop("Temperature in-base quantiles must contain 10th and 90th percentiles.\n")

  if("prec" %in% names(quantiles) && !all(c("q95", "q99") %in% names(quantiles$prec)))
    stop("Precipitation quantiles must contain 95th and 99th percentiles.\n")
}

get.temp.var.quantiles <- function(filled.data, date.series, bs.date.series, qtiles, bs.date.range, n, pad.data.with.first.last.values, in.base=FALSE) {
  base.data <- create.filled.series(filled.data, date.series, bs.date.series)
  if(in.base)
    return(list(outbase=zhang.running.qtile(base.data, dates.base=bs.date.series, qtiles=c(0.1, 0.9), bootstrap.range=bs.date.range, n=n, pad.data.with.first.last.values=pad.data.with.first.last.values),
                inbase=zhang.running.qtile(base.data, dates.base=bs.date.series, qtiles=c(0.1, 0.9), bootstrap.range=bs.date.range, n=n, pad.data.with.first.last.values=pad.data.with.first.last.values, get.bootstrap.data=TRUE)))
  else
    return(list(outbase=zhang.running.qtile(base.data, dates.base=bs.date.series, qtiles=c(0.1, 0.9), bootstrap.range=bs.date.range, n=n, pad.data.with.first.last.values=pad.data.with.first.last.values)))
}

get.prec.var.quantiles <- function(filled.prec, date.series, bs.date.range, qtiles=c(0.95, 0.99)) {
  wet.days <- !(is.na(filled.prec) | filled.prec < 1)
  inset <- date.series >= bs.date.range[1] & date.series <= bs.date.range[2] & !is.na(filled.prec) & wet.days
  pq <- quantile(filled.prec[inset], qtiles, type=8)
  names(pq) <- paste("q", qtiles * 100, sep="")
  return(pq)
}

#' Method for getting threshold quantiles for use in computing indices
#' 
#' This function creates threshold quantiles for use with climdexInput.raw
#' or climdexInput.csv.
#' 
#' This function takes input climate data at daily resolution, and produces as
#' output a set of threshold quantiles. This data structure can then be passed
#' to climdexInput.raw or climdexInput.csv.
#'
#' @param tmax Daily maximum temperature data.
#' @param tmin Daily minimum temperature data.
#' @param prec Daily total precipitation data.
#' @param tmax.dates Dates for the daily maximum temperature data.
#' @param tmin.dates Dates for the daily minimum temperature data.
#' @param prec.dates Dates for the daily total precipitation data.
#' @template climdexInput_common_params
#' @param quantiles Threshold quantiles for supplied variables.
#' @return A set of threshold quantiles
#' @note Units are assumed to be mm/day for precipitation and degrees Celsius
#' for temperature. No units conversion is performed internally.
#' 
#' @template climdexInput_raw_params_help
#' @template climdexInput_common_params_help
#' @examples
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
#' quantiles <- get.outofbase.quantiles(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#'
#' @export
get.outofbase.quantiles <- function(tmax=NULL, tmin=NULL, prec=NULL, tmax.dates=NULL, tmin.dates=NULL, prec.dates=NULL, base.range=c(1961, 1990), n=5, pad.data.with.first.last.values=FALSE, temp.qtiles=c(0.10, 0.90), prec.qtiles=c(0.95, 0.99)) {
  days.threshold <- 359
  check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n)
  
  all.dates <- c(tmin.dates, tmax.dates, prec.dates)
  last.day.of.year <- get.last.monthday.of.year(all.dates)
  cal <- attr(all.dates, "cal")

  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)
  new.date.range <- as.PCICt(paste(as.numeric(format(range(all.dates), "%Y", tz="GMT")), c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")

  bs.win.date.range <- get.bootstrap.windowed.range(bs.date.range, n)
  bs.date.series <- seq(bs.win.date.range[1], bs.win.date.range[2], by="day")

  quantiles <- list()

  if(!is.null(tmax)) {
    if(get.num.days.in.range(tmax.dates, bs.date.range) <= days.threshold)
      stop(paste("There is less than a year of tmax data within the base period. Consider revising your base range and/or check your input data."))
    filled.tmax <- create.filled.series(tmax, trunc(tmax.dates, "days"), date.series)
    quantiles$tmax <- get.temp.var.quantiles(filled.tmax, date.series, bs.date.series, temp.qtiles, bs.date.range, n, pad.data.with.first.last.values)
  } 

  if(!is.null(tmin)) {
    if(get.num.days.in.range(tmin.dates, bs.date.range) <= days.threshold)
      stop(paste("There is less than a year of tmin data within the base period. Consider revising your base range and/or check your input data."))
    filled.tmin <- create.filled.series(tmin, trunc(tmin.dates, "days"), date.series)
    quantiles$tmin <- get.temp.var.quantiles(filled.tmin, date.series, bs.date.series, temp.qtiles, bs.date.range, n, pad.data.with.first.last.values)
  }

  if(!is.null(prec)) {
    if(get.num.days.in.range(prec.dates, bs.date.range) <= days.threshold)
      stop(paste("There is less than a year of prec data within the base period. Consider revising your base range and/or check your input data."))
    filled.prec <- create.filled.series(prec, trunc(prec.dates, "days"), date.series)
    quantiles$prec <- get.prec.var.quantiles(filled.prec, date.series, bs.date.range, prec.qtiles)
  }
  return(quantiles)
}

#' Method for creating climdexInput object from vectors of data
#' 
#' This function creates a climdexInput object from data already ingested into
#' R.
#' 
#' This function takes input climate data at daily resolution, and produces as
#' output a ClimdexInput data structure. This data structure can then be passed
#' to any of the routines used to compute the Climdex indices. The indices
#' themselves are specified on the webpage cited in the references section.
#'
#' @template climdexInput_raw_help1 
#' @template climdexInput_raw_params
#' @template climdexInput_common_params
#' @param northern.hemisphere Whether this point is in the northern hemisphere.
#' @param quantiles Threshold quantiles for supplied variables.
#' @return An object of class \code{\link{climdexInput-class}} for use with
#' other climdex methods.
#' @note Units are assumed to be mm/day for precipitation and degrees Celsius
#' for temperature. No units conversion is performed internally.
#' 
#' @template climdexInput_raw_params_help
#' @template climdexInput_common_params_help
#' @examples
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
                             pad.data.with.first.last.values=FALSE, tavg=NULL, tavg.dates=NULL, quantiles=NULL, temp.qtiles=c(0.10, 0.90), prec.qtiles=c(0.95, 0.99)) {
  ## Make sure all of these arguments are valid...
  check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n, tavg, tavg.dates)

  all.dates <- c(tmin.dates, tmax.dates, prec.dates, tavg.dates)
  last.day.of.year <- get.last.monthday.of.year(all.dates)
  cal <- attr(all.dates, "cal")

  ## Convert base range (in years) to PCICt
  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)

  ## Get dates for normal data
  new.date.range <- as.PCICt(paste(as.numeric(format(range(all.dates), "%Y", tz="GMT")), c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")
  jdays <- get.jdays.replaced.feb29(get.jdays(date.series))
  
  ## Factors for dividing data up
  annual.factor <- factor(format(date.series, format="%Y", tz="GMT"))
  monthly.factor <- factor(format(date.series, format="%Y-%m", tz="GMT"))

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
  namask.ann <- lapply(filled.list, get.na.mask, annual.factor, 15)
  namask.mon <- lapply(filled.list, get.na.mask, monthly.factor, 3)
  names(namask.ann) <- names(namask.mon) <- names(filled.list)

  ## Pad data passed as base if we're missing endpoints...
  if(!have.quantiles) {
    quantiles <- environment()
    bs.win.date.range <- get.bootstrap.windowed.range(bs.date.range, n)
    bs.date.series <- seq(bs.win.date.range[1], bs.win.date.range[2], by="day")
    filled.list.base <- list(tmax=create.filled.series(filled.list$tmax, date.series, bs.date.series), tmin=create.filled.series(filled.list$tmin, date.series, bs.date.series))

    if(days.in.base['tmax'] > days.threshold)
      delayedAssign("tmax", get.temp.var.quantiles(filled.list$tmax, date.series, bs.date.series, temp.qtiles, bs.date.range, n, pad.data.with.first.last.values, TRUE), assign.env=quantiles)
    if(days.in.base['tmin'] > days.threshold)
      delayedAssign("tmin", get.temp.var.quantiles(filled.list$tmin, date.series, bs.date.series, temp.qtiles, bs.date.range, n, pad.data.with.first.last.values, TRUE), assign.env=quantiles)
    if(days.in.base['prec'] > days.threshold)
      delayedAssign("prec", get.prec.var.quantiles(filled.list$prec, date.series, bs.date.range, prec.qtiles), assign.env=quantiles)
  }
  
  return(new("climdexInput", data=filled.list, quantiles=quantiles, namask.ann=namask.ann, namask.mon=namask.mon, dates=date.series, jdays=jdays, base.range=bs.date.range, annual.factor=annual.factor, monthly.factor=monthly.factor, northern.hemisphere=northern.hemisphere))
}

#' Method for creating climdexInput object from CSV files
#' 
#' This function creates a climdexInput object from data in CSV files.
#' 
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
#' @return An object of class \code{\link{climdexInput-class}} for use with
#' other climdex methods.
#' @note Units are assumed to be mm/day for precipitation and degrees Celsius
#' for temperature. No units conversion is performed internally.
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
#' @template climdexInput_common_params_help
#' @examples
#' ## This would create a climdexInput object from a set of filenames (already
#' ## stored as variables), with a different date format.
#' \donttest{ci.csv <- climdexInput.csv(tmax.filename, tmin.filename,
#' prec.filename, date.types=list(list(fields=c("date"), format="%Y-%m-%d")))}
#'
#' @export
climdexInput.csv <- function(tmax.file=NULL, tmin.file=NULL, prec.file=NULL,
                             data.columns=list(tmin="tmin", tmax="tmax", prec="prec"), base.range=c(1961, 1990),
                             na.strings=NULL, cal="gregorian", date.types=NULL, n=5, northern.hemisphere=TRUE,
                             pad.data.with.first.last.values=FALSE, tavg.file=NULL, quantiles=NULL, temp.qtiles=c(0.10, 0.90), prec.qtiles=c(0.95, 0.99)) {
  get.and.check.data <- function(fn, datacol) {
    if(!is.null(fn)) {
      dat <- read.csv(fn, na.strings=na.strings)
      if(!(datacol %in% names(dat)))
        stop("Data column not found in tmin data.")
      return(list(dat=dat[!is.na(dat[,datacol]),],  dates=get.date.field(dat, cal, date.types)))
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
  
  return(climdexInput.raw(tmax=tmax$dat, tmin=tmin$dat, prec=prec$dat, tmax.dates=tmax$dates, tmin.dates=tmin$dates, prec.dates=prec$dates, base.range=base.range, n=n, northern.hemisphere=northern.hemisphere, pad.data.with.first.last.values=pad.data.with.first.last.values, tavg=tavg$dat, tavg.dates=tavg$dates, quantiles=quantiles, temp.qtiles=temp.qtiles, prec.qtiles=prec.qtiles))
}

#' Get date factors given frequency definition
#' 
#' This function takes a frequency definition and returns the appropriate date factor.
#' 
#' This is a convenience function which return the appropriate
#' date factor given an output time frequency (monthly or annual).
#' 
#' @param ci ClimdexInput object.
#' @param freq Frequency (either monthly or annual).
#' @return An NA mask or date factor for the given frequency.
#' @keywords ts climate
#' @examples
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
#' ## Get appropriate NA mask and date factor
#' date.factor.ann <- freq.to.factor(ci, "annual")
#' 
#' @export
freq.to.factor <- function(ci, freq) { switch(freq, annual=ci@annual.factor, monthly=ci@monthly.factor, NULL) }

#' Get NA mask given frequency definition
#' 
#' This function takes a frequency definition and returns the appropriate NA mask.
#' 
#' This is a convenience function which return the appropriate
#' NA mask given an output time frequency (monthly or annual).
#' 
#' @param ci ClimdexInput object.
#' @param freq Frequency (either monthly or annual).
#' @return An NA mask or date factor for the given frequency.
#' @keywords ts climate
#' @examples
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
#' ## Get appropriate NA mask and date factor
#' na.mask.ann <- freq.to.namask(ci, "annual")
#' 
#' @export
freq.to.namask <- function(ci, freq) { switch(freq, annual=ci@namask.ann, monthly=ci@namask.mon, NULL) }

#' Frost Days
#' 
#' This function computes the climdex index FD.
#' 
#' This function takes a climdexInput object as input and computes the FD (frost
#' days) climdex index: that is, the annual count of days where daily minimum
#' temperature drops below 0 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of frost days for each year.
#' @template generic_seealso_references
#' 
#' @templateVar cdxvar fd
#' @templateVar cdxdescription an annual timeseries of the number of frost days.
#' @template get_generic_example
#' 
#' @export
climdex.fd <- function(ci) { stopifnot(!is.null(ci@data$tmin)); return(number.days.op.threshold(ci@data$tmin, ci@annual.factor, 0, "<") * ci@namask.ann$tmin) }

#' Summer Days
#' 
#' This function computes the climdex index SU.
#' 
#' This function takes a climdexInput object as input and computes the SU (summer
#' days) climdex index: that is, the annual count of days where daily maximum
#' temperature exceeds 25 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of summer days for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar su
#' @templateVar cdxdescription an annual timeseries of the number of summer days.
#' @template get_generic_example
#' 
#' @export
climdex.su <- function(ci) { stopifnot(!is.null(ci@data$tmax)); return(number.days.op.threshold(ci@data$tmax, ci@annual.factor, 25, ">") * ci@namask.ann$tmax) }

#' Icing Days
#' 
#' This function computes the climdex index ID.
#' 
#' This function takes a climdexInput object as input and computes the ID (icing
#' days) climdex index: that is, the annual count of days where daily maximum
#' temperature is below 0 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of icing days for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar id
#' @templateVar cdxdescription an annual timeseries of the number of icing days.
#' @template get_generic_example
#' 
#' @export
climdex.id <- function(ci) { stopifnot(!is.null(ci@data$tmax)); return(number.days.op.threshold(ci@data$tmax, ci@annual.factor, 0, "<") * ci@namask.ann$tmax) }

#' Tropical Nights
#' 
#' This function computes the climdex index TR.
#' 
#' This function takes a climdexInput object as input and computes the TR
#' (tropical nights) climdex index: that is, the annual count of days where
#' daily minimum temperature stays above 20 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of frost days for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar tr
#' @templateVar cdxdescription an annual timeseries of the number of tropical nights.
#' @template get_generic_example
#' 
#' @export
climdex.tr <- function(ci) { stopifnot(!is.null(ci@data$tmin)); return(number.days.op.threshold(ci@data$tmin, ci@annual.factor, 20, ">") * ci@namask.ann$tmin) }

#' Growing Season Length
#' 
#' This function computes the growing season length (GSL) given the input.
#' 
#' This function takes a climdexInput object as input and computes the growing
#' season length based on this data.
#' 
#' Growing season length as defined by the climdex indices is the number of
#' days between the start of the first spell of warm days in the first half of
#' the year, and the start of the first spell of cold days in the second half
#' of the year. Spells of warm days are defined as six or more days with mean
#' temperature above 5 degrees Celsius; spells of cold days are defined as six
#' or more days with a mean temperature below 5 degrees Celsius.
#' 
#' The three alternate modes provided ('GSL_first', 'GSL_max', and 'GSL_sum')
#' are for testing purposes only. They differ considerably from the first
#' ('GSL') mode. All of them use a list of growing seasons -- here defined as
#' six or more consecutive days with a mean temperature greater than or equal
#' to 5 degrees Celsius, followed by either the end of the year or six or more
#' consecutive days with a mean temperature less than 5 degrees Celsius.
#' 'GSL_first' returns the first growing season found; 'GSL_max' returns the
#' longest growing season found; and 'GSL_sum' returns the total length of all
#' growing seasons found.
#' 
#' @param ci Object of type climdexInput.
#' @param gsl.mode Growing season length method to use.
#' @return A vector containing the number of days in the growing season for
#' each year.
#' @note Note that fclimdex results may differ from results using the first
#' ('GSL') mode due to bugs in fclimdex. Please ensure you are using the latest
#' version of fclimdex, as there have been numerous bug fixes and the results
#' should, at this point, match.
#' 
#' Please do not use the 'GSL_first', 'GSL_max', or 'GSL_sum' modes for
#' anything other than testing purposes at this time, nor should you rely on
#' this parameter being present in future versions of climdex.pcic.
#' @seealso \code{\link{growing.season.length}},
#' \code{\link{climdexInput.csv}}.
#' @references \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
#' @keywords ts climate
#' @templateVar cdxvar gsl
#' @templateVar cdxdescription an annual timeseries of the growing season length in days.
#' @template get_generic_example
#' 
#' @export
climdex.gsl <- function(ci, gsl.mode=c("GSL", "GSL_first", "GSL_max", "GSL_sum")) {
  stopifnot(!is.null(ci@data$tavg))
  ## Gotta shift dates so that July 1 is considered Jan 1 of same year in southern hemisphere
  if(ci@northern.hemisphere) {
    return(growing.season.length(ci@data$tavg, ci@annual.factor, ci@dates, ci@northern.hemisphere, gsl.mode=match.arg(gsl.mode)) * ci@namask.ann$tavg)
  } else {
    dates.POSIXlt <- as.POSIXlt(ci@dates)
    years <- dates.POSIXlt$year + 1900
    months <- dates.POSIXlt$mon + 1

    valid.years <- range(years)
    years.gsl <- years - floor((12 - months) / 6)

    inset <- years.gsl >= valid.years[1]
    gsl.factor <- factor(years.gsl[inset])
    gsl.temp.data <- ci@data$tavg[inset]
    namask.gsl <- get.na.mask(gsl.temp.data, gsl.factor, 15)
    namask.gsl[length(namask.gsl)] <- NA
    return((growing.season.length(gsl.temp.data, gsl.factor, ci@dates[inset], ci@northern.hemisphere, gsl.mode=match.arg(gsl.mode)) * namask.gsl))
  }
}

#' Monthly Maximum of Daily Maximum Temperature
#'
#' This function computes the climdex index TXx.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly or annual maximum of daily maximum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @return A vector containing the value of the index for each month.
#' @template generic_seealso_references
#' @templateVar cdxvar txx
#' @templateVar cdxdescription a monthly timeseries of maximum daily maximum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.txx <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax)); return(tapply.fast(ci@data$tmax, freq.to.factor(ci, match.arg(freq)), max) * freq.to.namask(ci, match.arg(freq))$tmax) }

#' Monthly Maximum of Daily Minimum Temperature
#'
#' This function computes the climdex index TNx.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly or annual maximum of daily minimum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @return A vector containing the value of the index for each month.
#' @template generic_seealso_references
#' @templateVar cdxvar tnx
#' @templateVar cdxdescription a monthly timeseries of maximum daily minimum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.tnx <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin)); return(tapply.fast(ci@data$tmin, freq.to.factor(ci, match.arg(freq)), max) * freq.to.namask(ci, match.arg(freq))$tmin) }

#' Monthly Minimum of Daily Maximum Temperature
#'
#' This function computes the climdex index TXn.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly or annual minimum of daily maximum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @return A vector containing the value of the index for each month.
#' @template generic_seealso_references
#' @templateVar cdxvar txn
#' @templateVar cdxdescription a monthly timeseries of minimum daily maximum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.txn <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax)); return(tapply.fast(ci@data$tmax, freq.to.factor(ci, match.arg(freq)), min) * freq.to.namask(ci, match.arg(freq))$tmax) }

#' Monthly Minimum of Daily Minimum Temperature
#'
#' This function computes the climdex index TNn.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly or annual minimum of daily minimum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @return A vector containing the value of the index for each month.
#' @template generic_seealso_references
#' @templateVar cdxvar tnn
#' @templateVar cdxdescription a monthly timeseries of minimum daily minimum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.tnn <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin)); return(tapply.fast(ci@data$tmin, freq.to.factor(ci, match.arg(freq)), min) * freq.to.namask(ci, match.arg(freq))$tmin) }

## Our implementation currently follows the example set by fclimdex for dealing with missing values, which is wrong; it biases results upwards when missing values are present.

#' Percent of Values Below 10th Percentile Daily Minimum Temperature
#' 
#' This function computes the climdex index TN10p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly or annual percent of values below the 10th percentile of baseline
#' daily minimum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tn10p
#' @templateVar cdxdescription a monthly timeseries of percentage of daily minimum temperature values which fall below the 10th percentile.
#' @template get_generic_example
#' 
#' @export
climdex.tn10p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(percent.days.op.threshold(ci@data$tmin, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmin$outbase$q10, ci@quantiles$tmin$inbase$q10, ci@base.range, "<") * freq.to.namask(ci, match.arg(freq))$tmin) }

#' Percent of Values Below 10th Percentile Daily Maximum Temperature
#' 
#' This function computes the climdex index TX10p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly or annual percent of values below the 10th percentile of baseline
#' daily maximum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tn10p
#' @templateVar cdxdescription a monthly timeseries of percentage of daily maximum temperature values which fall below the 10th percentile.
#' @template get_generic_example
#' 
#' @export
climdex.tx10p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmax$outbase$q10, ci@quantiles$tmax$inbase$q10, ci@base.range, "<") * freq.to.namask(ci, match.arg(freq))$tmax) }

#' Percent of Values Above 90th Percentile Daily Minimum Temperature
#' 
#' This function computes the climdex index TN90p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly or annual percent of values above the 90th percentile of baseline
#' daily minimum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tn10p
#' @templateVar cdxdescription a monthly timeseries of percentage of daily minimum temperature values which are above the 90th percentile.
#' @template get_generic_example
#' 
#' @export
climdex.tn90p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(percent.days.op.threshold(ci@data$tmin, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmin$outbase$q90, ci@quantiles$tmin$inbase$q90, ci@base.range, ">") * freq.to.namask(ci, match.arg(freq))$tmin) }

#' Percent of Values Above 90th Percentile Daily Maximum Temperature
#' 
#' This function computes the climdex index TX90p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly or annual percent of values above the 90th percentile of baseline
#' daily maximum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tn10p
#' @templateVar cdxdescription a monthly timeseries of percentage of daily maximum temperature values which are above the 90th percentile.
#' @template get_generic_example
#' 
#' @export
climdex.tx90p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmax$outbase$q90, ci@quantiles$tmax$inbase$q90, ci@base.range, ">") * freq.to.namask(ci, match.arg(freq))$tmax) }

#' Warm Spell Duration Index
#'
#' This function computes the climdex index WSDI.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index WSDI (Warm Spell Duration Index).
#' 
#' The warm spell duration index is defined as the number of days each year
#' which are part of a "warm spell". A "warm spell" is defined as a sequence of
#' 6 or more days in which the daily maximum temperature exceeds the 90th
#' percentile of daily maximum temperature for a 5-day running window
#' surrounding this day during the baseline period.
#' 
#' @template wcsdi_common
#' @seealso \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}},
#' \code{\link{threshold.exceedance.duration.index}}.
#' @references \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
#' @keywords ts climate
#' @templateVar cdxvar wsdi
#' @templateVar cdxdescription an annual timeseries of the warm spell duration index.
#' @template get_generic_example
#' 
#' @export
climdex.wsdi <- function(ci, spells.can.span.years=FALSE) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(threshold.exceedance.duration.index(ci@data$tmax, ci@annual.factor, ci@jdays, ci@quantiles$tmax$outbase$q90, ">", spells.can.span.years=spells.can.span.years) * ci@namask.ann$tmax) }

#' Cold Spell Duration Index
#' 
#' This function computes the climdex index CSDI.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index CSDI (Cold Spell Duration Index).
#'
#' The cold spell duration index is defined as the number of days
#' each year which are part of a "cold spell". A "cold spell" is defined as a
#' sequence of 6 or more days in which the daily minimum temperature is below
#' the 10th percentile of daily minimum temperature for a 5-day running window
#' surrounding this day during the baseline period.
#' 
#' @template wcsdi_common
#' @seealso \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}},
#' \code{\link{threshold.exceedance.duration.index}}.
#' @references \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
#' @keywords ts climate
#' @templateVar cdxvar csdi
#' @templateVar cdxdescription an annual timeseries of the cold spell duration index.
#' @template get_generic_example
#' 
#' @export
climdex.csdi <- function(ci, spells.can.span.years=FALSE) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(threshold.exceedance.duration.index(ci@data$tmin, ci@annual.factor, ci@jdays, ci@quantiles$tmin$outbase$q10, "<", spells.can.span.years=spells.can.span.years) * ci@namask.ann$tmax) }

#' Mean Diurnal Temperature Range
#' 
#' This function computes the diurnal temperature range on a monthly basis.
#' 
#' \code{climdex.dtr} computes the mean daily diurnal temperature range. The
#' frequency of observation can be either monthly or annual.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @return A vector containing the mean monthly or mean annual diurnal
#' temperature range.
#' @note This function creates results which may differ in the 3rd decimal
#' place from the results from fclimdex.
#' @template generic_seealso_references
#' @templateVar cdxvar dtr
#' @templateVar cdxdescription a monthly timeseries of mean diurnal temperature range.
#' @template get_generic_example
#' 
#' @export
climdex.dtr <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@data$tmax) && !is.null(ci@data$tavg)); return(mean.daily.temp.range(ci@data$tmax, ci@data$tmin, freq.to.factor(ci, match.arg(freq))) * freq.to.namask(ci, match.arg(freq))$tavg) }

#' Monthly Maximum 1-day Precipitation
#' 
#' This function computes the climdex index Rx1day.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index Rx1day: monthly or annual maximum 1-day precipitation.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @template rx5day_common
#' @template generic_seealso_references
#' @templateVar cdxvar rx1day
#' @templateVar cdxdescription a timeseries of monthly maximum 1-day precipitation.
#' @template get_generic_example
#' 
#' @export
climdex.rx1day <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$prec)); return(nday.consec.prec.max(ci@data$prec, freq.to.factor(ci, match.arg(freq)), 1) * freq.to.namask(ci, match.arg(freq))$prec) }

#' Monthly Maximum 5-day Consecutive Precipitation
#' 
#' This function computes the climdex index Rx5day.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index Rx5day: monthly or annual maximum 5-day consecutive precipitation.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param center.mean.on.last.day Whether to center the 5-day running mean on
#' the last day of the window, instead of the center day.
#' @template rx5day_common
#' @template generic_seealso_references
#' @templateVar cdxvar rx5day
#' @templateVar cdxdescription a timeseries of monthly maximum 5-day consecutive precipitation.
#' @template get_generic_example
#' 
#' @export
climdex.rx5day <- function(ci, freq=c("monthly", "annual"), center.mean.on.last.day=FALSE) { stopifnot(!is.null(ci@data$prec)); return(nday.consec.prec.max(ci@data$prec, freq.to.factor(ci, match.arg(freq)), 5, center.mean.on.last.day) * freq.to.namask(ci, match.arg(freq))$prec) }

#' Simple Precpitation Intensity Index
#' 
#' This function computes the climdex index SDII.
#' 
#' \code{climdex.sdii} computes the climdex index SDII, or Simple Precipitation
#' Intensity Index. This is defined as the sum of precipitation in wet days
#' (days with preciptitation over 1mm) during the year divided by the number of
#' wet days in the year.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the value of the index for each year.
#' @note fclimdex rounds to 1 decimal place, whereas climdex.sdii does not.
#' This results in some small differences.
#' @template generic_seealso_references
#' @templateVar cdxvar sdii
#' @templateVar cdxdescription a timeseries of annual SDII values.
#' @template get_generic_example
#' 
#' @export
climdex.sdii <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(simple.precipitation.intensity.index(ci@data$prec, ci@annual.factor) * ci@namask.ann$prec) }

#' Precipitation Exceeding 10mm Per Day
#' 
#' This function computes the climdex index R10mm.
#'
#' This function takes a climdexInput object as input and computes the climdex
#' index R10mm: the annual count of days where daily precipitation is more than 10mm per day.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the value of the index for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar r10mm
#' @templateVar cdxdescription an annual timeseries of the number of days where precipitation exceeds 10mm/day.
#' @template get_generic_example
#' 
#' @export
climdex.r10mm <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(number.days.op.threshold(ci@data$prec, ci@annual.factor, 10, ">=") * ci@namask.ann$prec) }

#' Precipitation Exceeding 20mm Per Day
#' 
#' This function computes the climdex index R20mm.
#'
#' This function takes a climdexInput object as input and computes the climdex
#' index R20mm: the annual count of days where daily precipitation is more than 20mm per day.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the value of the index for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar r20mm
#' @templateVar cdxdescription an annual timeseries of the number of days where precipitation exceeds 20mm/day.
#' @template get_generic_example
#' 
#' @export
climdex.r20mm <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(number.days.op.threshold(ci@data$prec, ci@annual.factor, 20, ">=") * ci@namask.ann$prec) }

#' Precipitation Exceeding A Specified Amount Per Day
#' 
#' This function computes the climdex index Rnnmm.
#'
#' This function takes a climdexInput object as input and computes the climdex
#' index Rnnmm: the annual count of days where daily precipitation is more than nn mm per day.
#' 
#' @param ci Object of type climdexInput.
#' @param threshold The threshold to be used for Rnnmm.
#' @return A vector containing the value of the index for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar rnnmm
#' @templateVar cdxdescription an annual timeseries of the number of days where precipitation exceeds 1 mm/day (the default).
#' @template get_generic_example
#' 
#' @export
climdex.rnnmm <- function(ci, threshold=1) {
  stopifnot(!is.null(ci@data$prec));
  if(!is.numeric(threshold) || length(threshold) != 1) stop("Please specify a single numeric threshold value.");

  return(number.days.op.threshold(ci@data$prec, ci@annual.factor, threshold, ">=") * ci@namask.ann$prec)
}

#' Maximum Consecutive Dry Days
#' 
#' This function computes the climdex index CDD.
#'
#' This function computes the climdex index CDD: the annual maximum length of dry spells, in days.
#' Dry spells are considered to be sequences of days where daily preciptation
#' is less than 1mm per day.
#' 
#' @template cdd_common
#' @templateVar cdxvar cdd
#' @templateVar cdxdescription an annual timeseries of the number of consecutive days where precipitation was less than 1mm/day.
#' @template get_generic_example
#' 
#' @export
climdex.cdd <- function(ci, spells.can.span.years=TRUE) { stopifnot(!is.null(ci@data$prec)); return(spell.length.max(ci@data$prec, ci@annual.factor, 1, "<", spells.can.span.years) * ci@namask.ann$prec) }

#' Maximum Consecutive Wet Days
#' 
#' This function computes the climdex index CWD.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index CWD: the annual maximum length of wet spells, in days.
#' Wet spells are considered to be sequences of days where daily precipitation
#' is at least 1mm per day.
#' 
#' @template cdd_common
#' @templateVar cdxvar cdd
#' @templateVar cdxdescription an annual timeseries of the number of consecutive days where precipitation was less than 1mm/day.
#' @template get_generic_example
#' 
#' @export
climdex.cwd <- function(ci, spells.can.span.years=TRUE) { stopifnot(!is.null(ci@data$prec)); return(spell.length.max(ci@data$prec, ci@annual.factor, 1, ">=", spells.can.span.years) * ci@namask.ann$prec) }

#' Total Daily Precipitation Exceeding 95\%ile Threshold
#' 
#' This function computes the climdex index R95pTOT.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index R95pTOT: the annual sum of precipitation in days where daily precipitation exceeds the
#' 95th percentile of daily precipitation in the base period.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing an annual timeseries of precipitation exceeding
#' the threshold.
#' @template generic_seealso_references
#' @templateVar cdxvar r95ptot
#' @templateVar cdxdescription an annual timeseries of the sum of precipitation where it exceeds the 95th percentile of precipitation in the base period.
#' @template get_generic_example
#' 
#' @export
climdex.r95ptot <- function(ci) { stopifnot(!is.null(ci@data$prec) && !is.null(ci@quantiles$prec)); return(total.precip.op.threshold(ci@data$prec, ci@annual.factor, ci@quantiles$prec['q95'], ">") * ci@namask.ann$prec) }

#' Total Daily Precipitation Exceeding 99\%ile Threshold
#' 
#' This function computes the climdex index R99pTOT.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index R99pTOT: the annual sum of precipitation in days where daily precipitation exceeds the
#' 99th percentile of daily precipitation in the base period.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing an annual timeseries of precipitation exceeding
#' the threshold.
#' @template generic_seealso_references
#' @templateVar cdxvar r99ptot
#' @templateVar cdxdescription an annual timeseries of the sum of precipitation where it exceeds the 99th percentile of precipitation in the base period.
#' @template get_generic_example
#' 
#' @export
climdex.r99ptot <- function(ci) { stopifnot(!is.null(ci@data$prec) && !is.null(ci@quantiles$prec)); return(total.precip.op.threshold(ci@data$prec, ci@annual.factor, ci@quantiles$prec['q99'], ">") * ci@namask.ann$prec) }

#' Total Daily Precipitation
#' 
#' This function computes the climdex index PRCPTOT.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index PRCPTOT: the annual sum of precipitation in wet days
#' (days where precipitation is at least 1mm).
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing an annual timeseries of precipitation in wet days.
#' @template generic_seealso_references
#' @templateVar cdxvar prcptot
#' @templateVar cdxdescription an annual timeseries of the sum of precipitation in wet days.
#' @template get_generic_example
#' 
#' @export
climdex.prcptot <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(total.precip.op.threshold(ci@data$prec, ci@annual.factor, 1, ">=") * ci@namask.ann$prec) }

all.indices <- c('fd', 'su', 'id', 'tr', 'gsl', 'txx', 'tnx', 'txn', 'tnn', 'tn10p', 'tx10p', 'tn90p', 'tx90p', 'wsdi', 'csdi',
                 'dtr', 'rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')

##
## HELPERS FINISHED. IMPLEMENTATION BELOW.
##

#' Get series length at ends
#' 
#' This function takes a series of boolean values and returns a list of
#' integers of the same length corresponding to the lengths at the ends of
#' sequences of TRUE values.
#' 
#' It can often be useful to know how long a series of boolean values is. This
#' function provides a method of knowing where and how long such sequences are.
#' 
#' @param x Sequence of booleans.
#' @param na.value Value to replace NAs with.
#' @return A vector consisting of the lengths of sequences of TRUE values at
#' the location of the last TRUE value in the sequence, and zeroes elsewhere.
#' @keywords ts climate
#' @examples
#' 
#' ## Get lengths of sequences of TRUE values in a sequence
#' series.lengths <- get.series.lengths.at.ends(c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE,
#' TRUE, TRUE, FALSE))
#' 
#' 
#' @export
get.series.lengths.at.ends <- function(x, na.value=FALSE) {
  n <- length(x)
  if(n == 1)
    return(as.numeric(x))

  res <- rep(0, n)
  x[is.na(x)] <- na.value

  ## Compare series to lag-1 and lag+1 series; false added to trigger state transition from TRUE at ends of series
  start <- which(x & !(c(FALSE, x[1:(n - 1)])))
  end <- which(x & !(c(x[2:n], FALSE)))
  res[end] <- end - start + 1
  return(res)
}

#' Number of days (less than, greater than, etc) a threshold
#' 
#' Produces sums of values that exceed (or are below) the specified threshold.
#' 
#' This function takes a data series, a threshold, an operator, and a factor to
#' aggregate by. It uses the operator to compare the threshold to the data
#' series, creating a series of booleans, then sums the booleans according to
#' the factor.
#' 
#' @param temp Sequence temperature values.
#' @param date.factor Factor to aggregate by.
#' @param threshold Threshold to use.
#' @param op Operator to use for comparison.
#' @return A vector consisting of the number of values that meet the criteria
#' in the given time period (as specified by \code{date.factor}).
#' @keywords ts climate
#' @examples
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
#' ## Calculate frost days.
#' fd <- number.days.op.threshold(ci@@data$tmin, ci@@annual.factor, 0, "<")
#' 
#' @export
number.days.op.threshold <- function(temp, date.factor, threshold, op="<") {
  stopifnot(is.numeric(c(temp, threshold)))
  return(tapply.fast(match.fun(op)(temp, threshold), date.factor, sum, na.rm=TRUE))
}

#' Flexible GSL function
#' 
#' This function computes the growing season length (GSL) given the input,
#' which is allowed to vary considerably from the ETCCDI definitions.
#' 
#' This function is the function used to implement \code{\link{climdex.gsl}}.
#' It's designed to be flexible to allow for experimentation and testing of new
#' thresholds and methods.
#' 
#' If you need to use this code for experimentation in the southern hemisphere,
#' you'll need to rip off the climdex.gsl code to rotate the year around so
#' that July 1st is treated as January 1st.
#' 
#' See \code{\link{climdex.gsl}} for more information on what \code{gsl.mode}
#' does.
#' 
#' @param daily.mean.temp Timeseries of daily mean temperature (in degrees C),
#' padded out to end on a year boundary (ie: starts on January 1st of some
#' year, ends on December 31st).
#' @param date.factor Factor of the same length as daily.mean.temp that divides
#' the timeseries up into years of data.
#' @param dates The corresponding series of dates.
#' @param northern.hemisphere Whether the data is from the northern hemisphere.
#' @param min.length The minimum number of days above or below the threshold
#' temperature that defines the start or end of a growing season.
#' @param t.thresh The temperature threshold for being considered part of a
#' growing season (in degrees C).
#' @param gsl.mode The growing season length mode (ETCCDI mode is "GSL").
#' @return A vector containing the number of days in the growing season for
#' each year.
#' @seealso \code{\link{climdex.gsl}}, \code{\link{climdexInput.csv}}.
#' @keywords ts climate
#' @examples
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
#' ## Create an annual timeseries of the growing season length in days.
#' gsl <- growing.season.length(ci@@data$tavg, ci@@annual.factor, ci@@dates,
#'                              ci@@northern.hemisphere, gsl.mode="GSL") * ci@@namask.ann$tavg
#' 
#' ## Print these out for testing purposes.
#' gsl
#' 
#' @export
growing.season.length <- function(daily.mean.temp, date.factor, dates, northern.hemisphere,
                                  min.length=6, t.thresh=5, gsl.mode=c("GSL", "GSL_first", "GSL_max", "GSL_sum")) {
  gsl.mode <- match.arg(gsl.mode)
  month.series <- get.months(dates)
  transition.month <- if(northern.hemisphere) 7 else 1
  if(gsl.mode == "GSL") {
    return(tapply.fast(1:length(daily.mean.temp), date.factor, function(idx) {
      temp.data <- daily.mean.temp[idx]
      ts.mid <- head(which(month.series[idx] == transition.month), n = 1)
      if(!length(ts.mid))
        return(NA)
      
      ts.len<- length(temp.data)
      gs.begin <- which(select.blocks.gt.length(temp.data[1:(ts.mid-1)] > t.thresh, min.length - 1))
      
      ## Growing season actually ends the day -before- the sequence of sketchy days
      gs.end <- which(select.blocks.gt.length(temp.data[ts.mid:ts.len] < t.thresh, min.length - 1)) - 1

      ## If no growing season start, 0 length; if no end, ends at end of year; otherwise, end - start + 1
      return(ifelse(length(gs.begin) == 0, 0, ifelse(length(gs.end) == 0, ts.len - gs.begin[1] + 1, gs.end[1] - gs.begin[1] + ts.mid)))
    }))
  } else {
    in.gsl <- !select.blocks.gt.length(!select.blocks.gt.length(daily.mean.temp >= t.thresh, min.length - 1), min.length - 1)
    warning("GSL_first, GSL_max, and GSL_sum are experimental alternative growing season length definitions. Use at your own risk.")
    
    innerfunc <- switch(gsl.mode, GSL_first=function(bl) { ifelse(any(bl > 0), (bl[bl > 0])[1], 0) }, GSL_max=max, GSL_sum=sum)
    return(tapply.fast(in.gsl, date.factor, function(ts) { block.lengths <- get.series.lengths.at.ends(ts); return(innerfunc(block.lengths)); }))
  }
}

## Sums up all of the arguments passed in
psum <- function (..., na.rm = FALSE) {
  args <- list(...)
  if(length(args) == 1 & is.list(args[[1]]))
    args <- args[[1]]
  stopifnot(length(unique(sapply(args, length))) == 1)
  if (na.rm)
    return(Reduce(function(x, y) { y[is.na(y)] <- 0; return(x + y); }, args, rep(0, length(args[[1]]))))
  else
    return(Reduce('+', args))
}

#' Lengths of strings of TRUE values
#' 
#' Computes fraction of days above or below the baseline threshold for each
#' day, and averages them using the date factor passed in.
#' 
#' This function computes fractions of days above or below baseline thresholds
#' for each day, then aggregates them using \code{date.factor}. It is used to
#' implement TN/TX 10/90p.
#' 
#' @param temp Sequence of temperature values.
#' @param dates Sequence of associated dates.
#' @param jdays Sequence of associated days of year.
#' @param date.factor Factor to aggregate data using.
#' @param threshold.outside.base Sequence of thresholds to be used for data
#' outside the base period.
#' @param base.thresholds Data structure containing sets of thresholds to be
#' used inside the base period; see \link{climdexInput-class}.
#' @param base.range Date range (type PCICt) of the baseline period.
#' @param op Comparison operator to use.
#' @return A vector consisting of the mean fraction of days above or below the
#' supplied set of thresholds.
#' @seealso \link{climdexInput-class}.
#' @keywords ts climate
#' @examples
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
#' ## Compute monthly tx90p.
#' tx90p <- percent.days.op.threshold(ci@@data$tmax, ci@@dates, ci@@jdays, ci@@monthly.factor,
#'                                    ci@@quantiles$tmax$outbase$q90,
#'                                    ci@@quantiles$tmax$inbase$q90,
#'                                    ci@@base.range, ">") * ci@@namask.mon$tmax
#' 
#' 
#' @export
percent.days.op.threshold <- function(temp, dates, jdays, date.factor, threshold.outside.base, base.thresholds, base.range, op='<') {
  f <- match.fun(op)
  dat <- f(temp, threshold.outside.base[jdays])
  
  inset <- dates >= base.range[1] & dates <= base.range[2]
  if(sum(inset) > 0) {
    jdays.base <- jdays[inset]
    years.base <- get.years(dates[inset])

    ## Get number of base years, subset temp data to base period only.
    temp.base <- temp[inset]
    years.base.range <- range(years.base)
    byrs <- (years.base.range[2] - years.base.range[1] + 1)

    ## Linearize thresholds, then compare them to the temperatures
    bdim <- dim(base.thresholds)
    dim(base.thresholds) <- c(bdim[1] * bdim[2], bdim[3])
    yday.byr.indices <- jdays.base + (years.base - get.years(base.range)[1]) * bdim[1]
    f.result <- f(rep(temp.base, byrs - 1), base.thresholds[yday.byr.indices,])
    dim(f.result) <- c(length(yday.byr.indices), bdim[3])

    ## Chop up data along the 2nd dim into a list; sum elements of the list
    dat[inset] <- psum(lapply(1:dim(f.result)[2], function(x) { f.result[,x] }), na.rm=TRUE) / (byrs - 1)
  }
  
  ret <- tapply.fast(dat, date.factor, function(x) { return(sum(x, na.rm=TRUE) / sum(!is.na(x))); }) * 100
  ret[is.nan(ret)] <- NA
  return(ret)
}

#' Sum of spell lengths exceeding daily threshold
#' 
#' This function returns the number of spells of more than \code{min.length}
#' days which exceed or are below the given threshold.
#' 
#' This routine compares data to the thresholds using the given operator,
#' generating a series of TRUE or FALSE values; these values are then filtered
#' to remove any sequences of less than \code{min.length} days of TRUE values.
#' It then computes the lengths of the remaining sequences of TRUE values
#' (spells) and sums their lengths.
#' 
#' The \code{spells.can.span.years} option controls whether spells must always
#' terminate at the end of a period, or whether they may continue until the
#' criteria ceases to be met or the end of the data is reached. The default for
#' fclimdex is FALSE.
#' 
#' @param daily.temp Data to compute index on.
#' @param date.factor Date factor to split by.
#' @param jdays Timeseries of days of year.
#' @param thresholds The thresholds to compare to.
#' @param op The operator to use to compare data to threshold.
#' @param min.length The minimum spell length to be considered.
#' @param spells.can.span.years Whether spells can span years.
#' @return A timeseries of maximum spell lengths for each period.
#' @seealso \code{\link{climdex.wsdi}}.
#' @keywords ts climate
#' @examples
#' 
#' prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
#' phony.date.factor <- factor(rep(1:2, each=5))
#' 
#' ## With spells spanning years...
#' alttedi <- threshold.exceedance.duration.index(prec.dat,
#' phony.date.factor, rep(1:5, 2), rep(1, 5), ">=", 2, TRUE)
#' 
#' ## Without spells spanning years...
#' tedi <- threshold.exceedance.duration.index(prec.dat, phony.date.factor,
#' rep(1:5, 2), rep(1, 5), ">=", 2, FALSE)
#' 
#' @export
threshold.exceedance.duration.index <- function(daily.temp, date.factor, jdays, thresholds, op=">", min.length=6, spells.can.span.years=TRUE) {
  stopifnot(is.numeric(c(daily.temp, thresholds, min.length)), is.factor(date.factor),
            is.function(match.fun(op)),
            min.length > 0)
  f <- match.fun(op)
  if(spells.can.span.years) {
    periods <- select.blocks.gt.length(f(daily.temp, thresholds[jdays]), min.length - 1)
    return(tapply.fast(periods, date.factor, sum))
  } else {
    ## fclimdex behaviour...
    return(tapply.fast(1:length(daily.temp), date.factor, function(idx) { sum(select.blocks.gt.length(f(daily.temp[idx], thresholds[jdays[idx]]), min.length - 1)) } ))
  }
}

## DTR
## Computes mean diurnal temperature range in each period (as specified by date.factor).
## Max and min temps are assumed to be same length
mean.daily.temp.range <- function(daily.max.temp, daily.min.temp, date.factor) {
  return(tapply.fast(daily.max.temp - daily.min.temp, date.factor, mean))
}

#' Number of days (less than, greater than, etc) a threshold
#' 
#' Produces sums of values that exceed (or are below) the specified threshold.
#' 
#' This function takes a data series, the number of days in the running window,
#' a date factor to aggregate by, and an optional modifier parameter
#' (center.mean.on.last.day). It computes the n-day running sum of
#' precipitation and returns the maximum n-day total precipitation per unit
#' time, as defined by \code{date.factor}.
#' 
#' @param daily.prec Daily timeseries of precipitation.
#' @param date.factor Factor to aggregate by.
#' @param ndays Number of days in the running window.
#' @param center.mean.on.last.day Whether to center the n-day running mean on
#' the last day of the series, instead of the middle day.
#' @return A vector consisting of the maximum n-day sum of precipitation per
#' time interval.
#' @keywords ts climate
#' @examples
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
#' ## Compute rx5day on a monthly basis.
#' rx5day <- nday.consec.prec.max(ci@@data$prec, ci@@monthly.factor, 5)
#' 
#' @export
nday.consec.prec.max <- function(daily.prec, date.factor, ndays, center.mean.on.last.day=FALSE) {
  if(ndays == 1) {
    return(tapply.fast(daily.prec, date.factor, max))
  } else {
    ## Ends of the data will be de-emphasized (padded with zero precip data); NAs replaced with 0
    new.series <- c(rep(0, floor(ndays / 2)), daily.prec, rep(0, floor(ndays / 2)))
    new.series[is.na(new.series)] <- 0
    prec.runsum <- runmean(new.series, k=ndays, endrule="trim") * ndays
    if(center.mean.on.last.day)
      prec.runsum <- c(rep(0, floor(ndays / 2)), prec.runsum[1:(length(prec.runsum) - floor(ndays / 2))])
    return(tapply.fast(prec.runsum, date.factor, max))
  }
}

#' Simple Precipitation Intensity Index
#' 
#' This function implements the ETCCDI Simple Precipitation Intensity Index.
#' 
#' The simple precipitation intensity index is computed by taking the sum of
#' precipitation in wet days (days with >1mm of precipitation), and dividing
#' that by the number of wet days in the period. This gives the mean
#' precipitation in wet days.
#' 
#' @param daily.prec Data to compute index on.
#' @param date.factor Date factor to split by.
#' @return The mean precipitation in wet days for each period (as defined by
#' date.factor).
#' @keywords ts climate
#' @examples
#' 
#' prec.dat <- c(0.1, 3.0, 4.3, 0.9, 1.3, 6.0, 0, 0, 4.0, 1)
#' phony.date.factor <- factor(rep(1:2, each=5))
#' sdii <- simple.precipitation.intensity.index(prec.dat, phony.date.factor)
#' 
#' @export
simple.precipitation.intensity.index <- function(daily.prec, date.factor) {
  return(tapply.fast(daily.prec, date.factor, function(prec) { idx <- prec >= 1 & !is.na(prec); if(sum(idx) == 0) { return(0); } else { return(sum(prec[idx], na.rm=TRUE) / sum(idx)) } } ))
}

#' Maximum spell length
#' 
#' This function returns the longest string of days which exceed or are below
#' the given threshold.
#' 
#' This routine compares data to the threshold using the given operator,
#' generating a series of TRUE or FALSE values. It then computes the lengths of
#' sequences of TRUE values (spells) and chooses the longest spell in each
#' period (as defined by date.factor).
#' 
#' The \code{spells.can.span.years} option controls whether spells must always
#' terminate at the end of a period, or whether they may continue until the
#' criteria ceases to be met or the end of the data is reached. The default for
#' fclimdex is TRUE.
#' 
#' @param daily.prec Data to compute index on.
#' @param date.factor Date factor to split by.
#' @param threshold The threshold to compare to.
#' @param op The operator to use to compare data to threshold.
#' @param spells.can.span.years Whether spells can span years.
#' @return A timeseries of maximum spell lengths for each period.
#' @seealso \code{\link{climdex.cdd}}.
#' @keywords ts climate
#' @examples
#' 
#' prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
#' phony.date.factor <- factor(rep(1:2, each=5))
#' 
#' ## With spells spanning years...
#' cwd <- spell.length.max(prec.dat, phony.date.factor, 1, ">=", TRUE)
#' 
#' ## Without spells spanning years...
#' altcwd <- spell.length.max(prec.dat, phony.date.factor, 1, ">=", FALSE)
#' 
#' @export
spell.length.max <- function(daily.prec, date.factor, threshold, op, spells.can.span.years) {
  bools <- match.fun(op)(daily.prec, threshold)

  if(spells.can.span.years) {
    all.true <- tapply.fast(bools, date.factor, all)
    max.spell <- tapply.fast(get.series.lengths.at.ends(bools), date.factor, max)
    
    ## Mask out values which are in the middle of a spell with NA
    na.mask <- c(1, NA)[as.integer((max.spell == 0) & all.true) + 1]
    return(max.spell * na.mask)
  } else {
    return(tapply.fast(bools, date.factor, function(x) { max(get.series.lengths.at.ends(x)) }))
  }
}

#' Sum of precipitation above a threshold
#' 
#' This function returns the sum of values above a threshold for each period
#' (as defined by date.factor).
#' 
#' This routine sums up all values which exceed or are below (depending on op)
#' the given threshold.
#' 
#' @param daily.prec Data to compute index on.
#' @param date.factor Date factor to split by.
#' @param threshold The threshold to compare to.
#' @param op The operator to use to compare data to threshold.
#' @return A timeseries of sums of numbers above the threshold for each period.
#' @seealso \code{\link{climdex.r99ptot}}.
#' @keywords ts climate
#' @examples
#' 
#' prec.dat <- c(0.1, 3.0, 4.3, 1.9, 1.3, 6.0, 0, 0, 4.0, 1)
#' phony.date.factor <- factor(rep(1:2, each=5))
#' 
#' ## Compute equiv of PRCPTOT
#' prec.sum <- total.precip.op.threshold(prec.dat, phony.date.factor, 1, ">=")
#' 
#' @export
total.precip.op.threshold <- function(daily.prec, date.factor, threshold, op) {
  f <- match.fun(op)
  return(tapply.fast(daily.prec, date.factor, function(pr) { return(sum(pr[f(pr, threshold)], na.rm=TRUE)) } ))
}

## Returns an n-day running quantile for each day of data (dimensions c(dpy, q))
## Data is assumed to be padded by floor(n/2) days on either end, and data is assumed to start on the (dpy - floor(n/2) + 1)'th day..
running.quantile <- function(data, n, q, dpy) {
  ret <- .Call("running_quantile_windowed", data, n, q, dpy, DUP=FALSE, PACKAGE='climdex.pcic')
  dim(ret) <- c(length(q), dpy)
  return(t(ret))
}

#' Select blocks of TRUE values of sufficient length.
#' 
#' Produces a sequence of booleans of the same length as input, with sequences
#' of TRUE values shorter than n replaced with FALSE.
#' 
#' This function takes a series of booleans and returns a sequence of booleans
#' of equal length, with all sequences of TRUE of length \code{n} or shorter
#' replaced with sequences of FALSE. NA values are replaced with
#' \code{na.value}.
#' 
#' @param d Sequence of booleans.
#' @param n Longest sequence of TRUE to replace with FALSE.
#' @param na.value Values to replace NAs with.
#' @return A vector of booleans, with the length \code{n} or less sequences of
#' TRUE replaced with FALSE.
#' @keywords ts climate
#' @examples
#' 
#' ## Return only the first sequence of TRUE... second sequence will be FALSE.
#' foo <- select.blocks.gt.length(c(rep(TRUE, 4), FALSE, rep(TRUE, 3)), 3)
#' 
#' @export
select.blocks.gt.length <- function(d, n, na.value=FALSE) {
  stopifnot(is.logical(d), is.numeric(n))

  if(n < 1)
    return(d)

  if(n >= length(d))
    return(rep(FALSE, length(d)))

  d[is.na(d)] <- na.value
  
  d2 <- Reduce(function(x, y) { return(c(rep(FALSE, y), d[1:(length(d) - y)]) & x) }, 1:n, d)
  return(Reduce(function(x, y) { return(c(d2[(y + 1):length(d2)], rep(FALSE, y)) | x) }, 1:n, d2))
}

#' Climdex quantile function
#' 
#' This function implements R's type=8 in a more efficient manner.
#' 
#' This is a reimplementation of R's type=8 created to improve the efficiency
#' of this package.
#' 
#' @param x Data to compute quantiles on.
#' @param q Quantiles to be computed.
#' @return A vector of the quantiles in question.
#' @seealso \code{\link{quantile}}
#' @keywords ts climate
#' @examples
#' 
#' ## Compute 10th, 50th, and 90th percentile of example data.
#' climdex.quantile(1:10, c(0.1, 0.5, 0.9))
#' 
#' @export
climdex.quantile <- function(x, q=c(0, 0.25, 0.5, 0.75, 1)) {
  return(.Call("c_quantile2", as.double(x), q, PACKAGE='climdex.pcic'))
}
