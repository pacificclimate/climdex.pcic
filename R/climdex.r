library(caTools)
library(PCICt)

bs.pctile.names <- c("tx10thresh", "tx90thresh", "tn10thresh", "tn90thresh")
prec.pctile.names <- c("r95thresh", "r99thresh")

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
  errors <- c()

  if(!all(sapply(c("tmax", "tmin", "tavg", "prec", "dates", "jdays", "annual.factor", "monthly.factor"), function(y) { length(slot(x, y)) }) == length(x@tmax)))
    errors <- c(errors, "Data fields, dates, and factors must all be of the same length")

  ## Check that namask.mon and namask.ann have columns for each of the variables
  if(sum(c("tmax", "tmin", "tavg", "prec") %in% names(x@namask.ann)) != 4)
    errors <- c(errors, "NA mask for annual must contain data for tmax, tmin, tavg, and prec.")

  if(sum(c("tmax", "tmin", "tavg", "prec") %in% names(x@namask.mon)) != 4)
    errors <- c(errors, "NA mask for monthly must contain data for tmax, tmin, tavg, and prec.")

  ## Check that running.pctile.base running.pctile.notbase, and pctile all contain the proper named data.
  for(n in bs.pctile.names)
    if(!(n %in% names(x@running.pctile.base)) || (is.null(x@running.pctile.base[[n]]) && get.num.days.in.range(x@dates, x@base.range) > 0))
      errors <- c(errors, paste("running.pctile.base$", n, " is null but the data contains dates in the base period; or it does not contain ", n, " at all.", sep=""))

  if(sum(bs.pctile.names %in% names(x@running.pctile.notbase)) != length(bs.pctile.names))
    errors <- c(errors, "running.pctile.notbase does not contain at least one of tx10thresh, tx90thresh, tn10thresh, and tn90thresh.")

  if(sum(prec.pctile.names %in% names(x@prec.pctile)) != length(prec.pctile.names))
    errors <- c(errors, "pctile does not contain at least one of r95pthresh and r99pthresh.")

  if(length(x@northern.hemisphere) != 1)
    errors <- c(errors, "northern.hemisphere must be of length 1.")
  
  if(length(errors) == 0)
    return(TRUE)
  else
    return(errors)
}

## Class definition declaration
setClass("climdexInput",
         representation(tmax = "numeric",
                        tmin = "numeric",
                        tavg = "numeric",
                        prec = "numeric",
                        namask.ann = "data.frame",
                        namask.mon = "data.frame",
                        running.pctile.base = "list",
                        running.pctile.notbase = "data.frame",
                        prec.pctile = "numeric",
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
  as.integer(format(dates, "%j", tz="GMT"))
}

## Get year
get.years <- function(dates) {
  as.integer(format(dates, format="%Y", tz="GMT"))
}

## Get month number
get.months <- function(dates) {
  as.integer(format(dates, format="%m", tz="GMT"))
}

## Juggle the list so that day 366 == day 365
get.jdays.replaced.feb29 <- function(dates) {
  jdays <- get.jdays(dates)
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

## Do the Zhang boostrapping method described in Xuebin Zhang et al's 2005 paper, "Avoiding Inhomogeneity in Percentile-Based Indices of Temperature Extremes" J.Clim vol 18 pp.1647-1648, "Removing the 'jump'".
## Expects PCICt for all dates
zhang.bootstrap.qtile <- function(x, dates, qtiles, bootstrap.range, include.mask=NULL, n=5, pad.data.with.first.last.values=FALSE) {
  window <- floor(n / 2)

  dpy <- ifelse(is.null(attr(dates, "dpy")), 365, attr(dates, "dpy"))
  inset <- get.bootstrap.set(dates, bootstrap.range, n)
  nyears <- floor(sum(inset) / dpy)
  
  if(!is.null(include.mask))
    x[include.mask] <- NA

  bs.data <- x[inset]
  if(pad.data.with.first.last.values) {
    bs.data[1:window] <- bs.data[window + 1]
    bs.data[length(bs.data) - 0:(window - 1)] <- bs.data[length(bs.data) - window]
  }
  
  ## This routine is written as described in Zhang et al, 2005 as referenced above.
  yidx <- 1:nyears
  d <- sapply(yidx, function(idx.to.omit) {
    ## Index is computed assuming `window` leading days of data...
    omit.index <- window + (1:dpy) + ((idx.to.omit - 1) * dpy)
    return(sapply(yidx[yidx != idx.to.omit], function(idx.to.replace.with) {
      replace.index <- window + (1:dpy) + ((idx.to.replace.with - 1) * dpy)
      bs.data[omit.index] <- bs.data[replace.index]
      return(running.quantile(bs.data, n, qtiles, dpy))
    }, simplify="array"))
  }, simplify="array" )

  d <- aperm(d, perm=c(1, 4, 3, 2))
  ## new dims: dpy, nyears, nyears-1, length(quantiles)
  
  return(lapply(1:length(qtiles), function(x) { d[,,,x] }))
}

## Calculate a running quantile on the data set over the bootstrap range
zhang.running.qtile <- function(x, dates, dates.base, qtiles, bootstrap.range, include.mask=NULL, n=5, pad.data.with.first.last.values=FALSE) {
  inset <- get.bootstrap.set(dates.base, bootstrap.range, n)
  dpy <- ifelse(is.null(attr(dates, "dpy")), 365, attr(dates, "dpy"))
  
  if(!is.null(include.mask))
    x[include.mask] <- NA

  bs.data <- x[inset]
  window <- floor(n / 2)
  if(pad.data.with.first.last.values) {
    bs.data[1:window] <- bs.data[window + 1]
    bs.data[length(bs.data) - 0:(window - 1)] <- bs.data[length(bs.data) - window]
  }

  running.quantile(bs.data, n, qtiles, dpy)
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
  if(!(inherits(tmax.dates, "PCICt") && inherits(tmin.dates, "PCICt") && inherits(prec.dates, "PCICt") && (is.null(tavg.dates) || inherits(tavg.dates, "PCICt"))))
    stop("Dates must all be of class PCICt.")

  ## FIXME: Add type checks on all args 

  if(is.null(tavg) != is.null(tavg.dates))
    stop("If passing in tavg, must pass in tavg dates too..")

  if(length(tmin) != length(tmin.dates))
    stop("Length of tmin data and tmin dates do not match.")
  
  if(length(tmax) != length(tmax.dates))
    stop("Length of tmax data and tmax dates do not match.")
  
  if(length(prec) != length(prec.dates))
    stop("Length of prec data and prec dates do not match.")

  if(!is.null(tavg) && length(tavg) != length(tavg.dates)) 
    stop("Length of tavg data and tavg dates do not match.")

  if(!(length(base.range) == 2 && is.numeric(base.range)))
    stop("Invalid base date range; expecting vector of 2 numeric years.")

  dates.list <- list(tmax.dates, tmin.dates, prec.dates)
  if(!is.null(tavg))
    dates.list <- c(dates.list, tavg.dates)
  if(any(!sapply(dates.list, inherits, "PCICt")))
    stop("Dates must be of class PCICt.")

  if(n != 5)
    warning("Use of n != 5 varies from the Climdex definition. Use at your own risk.")
}

## Check validity of quantile input.
check.quantile.validity <- function(temp.quantiles.notbase, temp.quantiles.base, prec.quantiles, all.dates, base.range) {
  ## Check that all is well with out-of-base quantiles
  if(is.null(temp.quantiles.notbase) != is.null(prec.quantiles))
    stop("Either all of the out-of-base quantiles must be passed in, or none.")
  
  ## Check that tmin/tmax percentiles (base and not) contain 10th and 90th with names as given
  if(!is.null(temp.quantiles.notbase) && !all(bs.pctile.names %in% names(temp.quantiles.notbase)))
    stop("temp.quantiles.notbase must contain tx10thresh (10th percentile), tx90thresh (90th percentile), tn10thresh (10th percentile) and tn90thresh (90th percentile).")
  ## FIXME: Add dimensionality checks
  
  ## Check that prec quantiles contain 95th and 99th
  if(!is.null(prec.quantiles) && !all(prec.pctile.names %in% names(prec.quantiles)))
    stop("Precipitation quantiles do not contain r95thresh (95th percentile of wet days) and r99thresh (99th percentile of wet days).")

  if(!is.null(temp.quantiles.base)) {
    if(!all(bs.pctile.names %in% names(temp.quantiles.base)))
      stop("temp.quantiles.base must contain tx10thresh (10th percentile), tx90thresh (90th percentile), tn10thresh (10th percentile) and tn90thresh (90th percentile).")
    ## FIXME: Add dimensionality checks
  }
}

## Get quantiles for out-of-base period; for use when computing indices on future data using historical quantiles.
get.outofbase.quantiles <- function(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range=c(1961, 1990), n=5, pad.data.with.first.last.values=FALSE) {
  check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n)

  ## Get last day of year (yay variable calendars)
  cal <- attr(tmax.dates, "cal")
  last.day.of.year <- "12-31"
  if(!is.null(attr(tmax.dates, "months")))
    last.day.of.year <- paste("12", attr(tmax.dates, "months")[12], sep="-")

  ## Convert base range (in years) to PCICt
  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)
  
  ## Get dates for normal data
  all.dates <- c(tmin.dates, tmax.dates, prec.dates)
  new.date.range <- as.PCICt(paste(as.numeric(format(range(all.dates), "%Y", tz="GMT")), c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")

  ## Filled data...
  filled.tmax <- create.filled.series(tmax, trunc(tmax.dates, "days"), date.series)
  filled.tmin <- create.filled.series(tmin, trunc(tmin.dates, "days"), date.series)
  filled.prec <- create.filled.series(prec, trunc(prec.dates, "days"), date.series)
  
  ## Establish some truth values for later use in logic...
  days.threshold <- 359
  if(!all(sapply(list(tmax.dates, tmin.dates, prec.dates), get.num.days.in.range, bs.date.range) > days.threshold))
    stop("There is less than a year of data for at least one of tmin, tmax, and prec. Consider revising your base range and/or check your input data.")

  ## DeMorgan's laws FTW
  wet.days <- !(is.na(filled.prec) | filled.prec < 1)

  ## Pad data passed as base if we're missing endpoints...
  bs.win.date.range <- get.bootstrap.windowed.range(bs.date.range, n)
  bs.date.series <- seq(bs.win.date.range[1], bs.win.date.range[2], by="day")
  filled.list.base <- list(tmax=create.filled.series(filled.tmax, date.series, bs.date.series), tmin=create.filled.series(filled.tmin, date.series, bs.date.series))
  
  bs.pctile <- do.call(data.frame, lapply(filled.list.base[1:2], zhang.running.qtile, dates=date.series, dates.base=bs.date.series, qtiles=c(0.1, 0.9), bootstrap.range=bs.date.range, n=n, pad.data.with.first.last.values=pad.data.with.first.last.values))
  inset <- date.series >= bs.date.range[1] & date.series <= bs.date.range[2] & !is.na(filled.prec) & wet.days
  prec.pctile <- quantile(filled.prec[inset], c(0.95, 0.99), type=8)
  names(bs.pctile) <- c("tx10thresh", "tx90thresh", "tn10thresh", "tn90thresh")
  names(prec.pctile) <- c("r95thresh", "r99thresh")
  return(list(bs.pctile=bs.pctile, prec.pctile=prec.pctile))
}

## Creates climdexInput data structure given input.
## Temp. units: degrees C, Prec. units: mm per day
climdexInput.raw <- function(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range=c(1961, 1990), n=5,
                             northern.hemisphere=TRUE
                             , pad.data.with.first.last.values=FALSE,
                             tavg=NULL, tavg.dates=NULL, temp.quantiles.notbase=NULL,
                             prec.quantiles=NULL, temp.quantiles.base=NULL) {
  ## Make sure all of these arguments are valid...
  check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n, tavg, tavg.dates)

  all.dates <- c(tmin.dates, tmax.dates, prec.dates)
  if(!is.null(tavg.dates))
    all.dates <- c(all.dates, tavg.dates)

  ## Get last day of year (yay variable calendars)
  cal <- attr(tmax.dates, "cal")
  last.day.of.year <- "12-31"
  if(!is.null(attr(tmax.dates, "months")))
    last.day.of.year <- paste("12", attr(tmax.dates, "months")[12], sep="-")

  ## Convert base range (in years) to PCICt
  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)

  ## Check that quantiles are valid
  check.quantile.validity(temp.quantiles.notbase, temp.quantiles.base, prec.quantiles, all.dates, bs.date.range)
  
  ## Get dates for normal data
  new.date.range <- as.PCICt(paste(as.numeric(format(range(all.dates), "%Y", tz="GMT")), c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")
  jdays <- get.jdays.replaced.feb29(date.series)
  
  ## Factors for dividing data up
  annual.factor <- as.factor(format(date.series, "%Y", tz="GMT"))
  monthly.factor <- as.factor(format(date.series, "%Y-%m", tz="GMT"))

  ## Filled data...
  filled.tmax <- create.filled.series(tmax, trunc(tmax.dates, "days"), date.series)
  filled.tmin <- create.filled.series(tmin, trunc(tmin.dates, "days"), date.series)
  filled.prec <- create.filled.series(prec, trunc(prec.dates, "days"), date.series)
  filled.tavg <- (filled.tmax + filled.tmin) / 2
  if(!is.null(tavg))
    filled.tavg <- create.filled.series(tavg, trunc(tavg.dates, "days"), date.series)
  
  ## Establish some truth values for later use in logic...
  days.threshold <- 359
  data.in.base.period <- all(sapply(list(tmax.dates, tmin.dates, prec.dates), get.num.days.in.range, bs.date.range) != 0)
  have.notbase.quantiles <- !is.null(temp.quantiles.notbase)
  have.base.quantiles <- !is.null(temp.quantiles.base)
  
  if(!all(sapply(list(tmax.dates, tmin.dates, prec.dates), get.num.days.in.range, bs.date.range) > days.threshold)) {
    if(any(sapply(list(tmax.dates, tmin.dates, prec.dates), get.num.days.in.range, bs.date.range) > 0) && !have.base.quantiles)
      stop("Base quantiles were not specified and there is insufficient data to compute them. Consider revising your base range and/or check your input data.")
    if(!have.notbase.quantiles) {
      stop("Out of base quantiles were not specified and there is insufficient data to compute them. Consider revising your base range and/or check your input data.")
    }
  }

  filled.list <- list(tmax=filled.tmax, tmin=filled.tmin, tavg=filled.tavg, prec=filled.prec)
  filled.list.names <- names(filled.list)

  ## NA masks
  namask.ann <- do.call(data.frame, lapply(filled.list, get.na.mask, annual.factor, 15))
  namask.mon <- do.call(data.frame, lapply(filled.list, get.na.mask, monthly.factor, 3))
  colnames(namask.ann) <- colnames(namask.mon) <- filled.list.names

  ## DeMorgan's laws FTW
  wet.days <- !(is.na(filled.prec) | filled.prec < 1)

  ## Pad data passed as base if we're missing endpoints...
  bs.pctile.base <- temp.quantiles.base
  bs.pctile <- temp.quantiles.notbase
  prec.pctile <- prec.quantiles
  if((!have.base.quantiles && data.in.base.period) || !have.notbase.quantiles) {
    bs.win.date.range <- get.bootstrap.windowed.range(bs.date.range, n)
    bs.date.series <- seq(bs.win.date.range[1], bs.win.date.range[2], by="day")
    filled.list.base <- list(tmax=create.filled.series(filled.tmax, date.series, bs.date.series), tmin=create.filled.series(filled.tmin, date.series, bs.date.series))

    if(!have.notbase.quantiles) {
      bs.pctile <- do.call(data.frame, lapply(filled.list.base[1:2], zhang.running.qtile, dates=date.series, dates.base=bs.date.series, qtiles=c(0.1, 0.9), bootstrap.range=bs.date.range, n=n, pad.data.with.first.last.values=pad.data.with.first.last.values))
      inset <- date.series >= bs.date.range[1] & date.series <= bs.date.range[2] & !is.na(filled.prec) & wet.days
      prec.pctile <- quantile(filled.prec[inset], c(0.95, 0.99), type=8)
      names(bs.pctile) <- bs.pctile.names
      names(prec.pctile) <- prec.pctile.names
    }

    if(!have.base.quantiles) {
      bs.pctile.base <- do.call(c, lapply(filled.list.base[1:2], zhang.bootstrap.qtile, bs.date.series, c(0.1, 0.9), bs.date.range, n=n, pad.data.with.first.last.values=pad.data.with.first.last.values))
      names(bs.pctile.base) <- bs.pctile.names
    }
  }

  ## Create a sane list here which won't blow things up when bits of bs.pctile.base are passed as args to various helper funcs.
  if(is.null(bs.pctile.base))
    bs.pctile.base <- list(tx10thresh=NULL, tx90thresh=NULL, tn10thresh=NULL, tn90thresh=NULL)
  
  return(new("climdexInput", tmax=filled.tmax, tmin=filled.tmin, tavg=filled.tavg, prec=filled.prec, namask.ann=namask.ann, namask.mon=namask.mon, running.pctile.base=bs.pctile.base, running.pctile.notbase=bs.pctile, prec.pctile=prec.pctile, dates=date.series, jdays=jdays, base.range=bs.date.range, annual.factor=annual.factor, monthly.factor=monthly.factor, northern.hemisphere=northern.hemisphere))
}

## Creates climdexInput data structure given input CSV files consisting of date columns and data.
## Temp. units: degrees C, Prec. units: mm per day
climdexInput.csv <- function(tmax.file, tmin.file, prec.file,
                             data.columns=list(tmin="tmin", tmax="tmax", prec="prec"), base.range=c(1961, 1990),
                             na.strings=NULL, cal="gregorian", date.types=NULL, n=5, northern.hemisphere=TRUE,
                             pad.data.with.first.last.values=FALSE, tavg.file=NULL, temp.quantiles.notbase=NULL,
                             prec.quantiles=NULL, temp.quantiles.base=NULL) {
  
  if(sum(c("tmax", "tmin", "prec") %in% names(data.columns)) != 3)
    stop("Must specify names of all data columns (tmin, tmax, prec).")

  if(!is.null(tavg.file) && !("tavg" %in% names(data.columns)))
    stop("If specifying tavg, must also specify the name of the column tavg data is in.")
  
  if(missing(date.types))
    date.types <- list(list(fields=c("year", "jday"), format="%Y %j"),
                       list(fields=c("year", "month", "day"), format="%Y %m %d"))
  else
    if(any(!sapply(date.types, function(x) { return(sum(c("fields", "format") %in% names(x)) == 2 && is.character(x$fields) && is.character(x$format)) } )))
      stop("Invalid date.types specified. See ?climdexInput.csv .")
  
  tmin.dat <- read.csv(tmin.file, na.strings=na.strings)
  tmax.dat <- read.csv(tmax.file, na.strings=na.strings)
  prec.dat <- read.csv(prec.file, na.strings=na.strings)

  ## Optionally read in tavg
  tavg.dat <- NULL
  if(!is.null(tavg.file))
    tavg.dat <- read.csv(tavg.file, na.strings=na.strings)

  ## Check for data columns in data files
  if(!(data.columns$tmin %in% names(tmin.dat)))
    stop("Data column not found in tmin data.")
  if(!(data.columns$tmax %in% names(tmax.dat)))
    stop("Data column not found in tmax data.")
  if(!(data.columns$prec %in% names(prec.dat)))
    stop("Data column not found in prec data.")
  if(!is.null(tavg.dat) & !(data.columns$tavg %in% names(tavg.dat)))
    stop("Data column not found in tavg data.")

  ## This is to deal with fclimdex's broken input data, which includes February 31st, amongst other joyous things
  tmin.dat <- tmin.dat[!is.na(tmin.dat[,data.columns$tmin]),]
  tmax.dat <- tmax.dat[!is.na(tmax.dat[,data.columns$tmax]),]
  prec.dat <- prec.dat[!is.na(prec.dat[,data.columns$prec]),]
  
  tmin.dates <- get.date.field(tmin.dat, cal, date.types)
  tmax.dates <- get.date.field(tmax.dat, cal, date.types)
  prec.dates <- get.date.field(prec.dat, cal, date.types)

  ## Handle tavg being optionally present -- make sure we don't have to do anything special later.
  if(!is.null(tavg.dat)) {
    tavg.dat <- tavg.dat[!is.na(tavg.dat[,data.columns$tavg]),]
    tavg.dates <- get.date.field(tavg.dat, cal, date.types)
    tavg.input <- tavg.dat[,data.columns$tavg]
  } else {
    tavg.dates <- NULL
    tavg.input <- NULL
  }
  
  return(climdexInput.raw(tmax=tmax.dat[,data.columns$tmax], tmin=tmin.dat[,data.columns$tmin], prec=prec.dat[,data.columns$prec], tmax.dates=tmax.dates, tmin.dates=tmin.dates, prec.dates=prec.dates, base.range=base.range, n=n, northern.hemisphere=northern.hemisphere, pad.data.with.first.last.values=pad.data.with.first.last.values, tavg=tavg.input, tavg.dates=tavg.dates, temp.quantiles.notbase=temp.quantiles.notbase, prec.quantiles=prec.quantiles, temp.quantiles.base=temp.quantiles.base))
}

## Returns appropriate factor for given time frequency
freq.to.factor <- function(ci, freq) { switch(freq, annual=ci@annual.factor, monthly=ci@monthly.factor, NULL) }

## Returns appropriate NA mask for given time frequency
freq.to.namask <- function(ci, freq) { switch(freq, annual=ci@namask.ann, monthly=ci@namask.mon, NULL) }

## FD: Differences of 1-4 days in some years.
climdex.fd <- function(ci) { return(number.days.op.threshold(ci@tmin, ci@annual.factor, 0, "<") * ci@namask.ann$tmin) }

## SU: Differences of 1-4 days in some years.
climdex.su <- function(ci) { return(number.days.op.threshold(ci@tmax, ci@annual.factor, 25, ">") * ci@namask.ann$tmax) }

## ID: Differences of 1-4 days in some years.
climdex.id <- function(ci) { return(number.days.op.threshold(ci@tmax, ci@annual.factor, 0, "<") * ci@namask.ann$tmax) }

## TR: Differences of 1-4 days in some years.
climdex.tr <- function(ci) { return(number.days.op.threshold(ci@tmin, ci@annual.factor, 20, ">") * ci@namask.ann$tmin) }

## GSL: Modes other than "GSL" may not be well tested.
climdex.gsl <- function(ci, gsl.mode=c("GSL", "GSL_first", "GSL_max", "GSL_sum")) {
  ## Gotta shift dates so that July 1 is considered Jan 1 of same year in southern hemisphere
  if(ci@northern.hemisphere) {
    return(growing.season.length(ci@tavg, ci@annual.factor, ci@dates, ci@northern.hemisphere, gsl.mode=match.arg(gsl.mode)) * ci@namask.ann$tavg)
  } else {
    valid.date.range <- range(ci@dates)
    valid.years <- get.years(valid.date.range)
    years.gsl <- get.years(ci@dates) - floor((12 - get.months(ci@dates)) / 6)

    inset <- years.gsl >= valid.years[1]
    gsl.factor <- factor(as.character(years.gsl[inset]))
    gsl.temp.data <- ci@tavg[inset]
    namask.gsl <- get.na.mask(gsl.temp.data, gsl.factor, 15)
    namask.gsl[length(namask.gsl)] <- NA
    return((growing.season.length(gsl.temp.data, gsl.factor, ci@dates[inset], ci@northern.hemisphere, gsl.mode=match.arg(gsl.mode)) * namask.gsl))
  }
}

## TXx
climdex.txx <- function(ci, freq=c("monthly", "annual")) { return(tapply.fast(ci@tmax, freq.to.factor(ci, match.arg(freq)), max) * freq.to.namask(ci, match.arg(freq))$tmax) }

## TNx
climdex.tnx <- function(ci, freq=c("monthly", "annual")) { return(tapply.fast(ci@tmin, freq.to.factor(ci, match.arg(freq)), max) * freq.to.namask(ci, match.arg(freq))$tmin) }

## TXn
climdex.txn <- function(ci, freq=c("monthly", "annual")) { return(tapply.fast(ci@tmax, freq.to.factor(ci, match.arg(freq)), min) * freq.to.namask(ci, match.arg(freq))$tmax) }

## TNn
climdex.tnn <- function(ci, freq=c("monthly", "annual")) { return(tapply.fast(ci@tmin, freq.to.factor(ci, match.arg(freq)), min) * freq.to.namask(ci, match.arg(freq))$tmin) }

## TN10p
## Our implementation currently follows the example set by fclimdex for dealing with missing values, which is wrong; it biases results upwards when missing values are present.
climdex.tn10p <- function(ci, freq=c("monthly", "annual")) { return(percent.days.op.threshold(ci@tmin, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@running.pctile.notbase$tn10thresh, ci@running.pctile.base$tn10thresh, ci@base.range, "<") * freq.to.namask(ci, match.arg(freq))$tmin) }

## TX10p
climdex.tx10p <- function(ci, freq=c("monthly", "annual")) { return(percent.days.op.threshold(ci@tmax, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@running.pctile.notbase$tx10thresh, ci@running.pctile.base$tx10thresh, ci@base.range, "<") * freq.to.namask(ci, match.arg(freq))$tmax) }

## TN90p
climdex.tn90p <- function(ci, freq=c("monthly", "annual")) { return(percent.days.op.threshold(ci@tmin, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@running.pctile.notbase$tn90thresh, ci@running.pctile.base$tn90thresh, ci@base.range, ">") * freq.to.namask(ci, match.arg(freq))$tmin) }

## TX90p
climdex.tx90p <- function(ci, freq=c("monthly", "annual")) { return(percent.days.op.threshold(ci@tmax, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@running.pctile.notbase$tx90thresh, ci@running.pctile.base$tx90thresh, ci@base.range, ">") * freq.to.namask(ci, match.arg(freq))$tmax) }

## WSDI
climdex.wsdi <- function(ci, spells.can.span.years=FALSE) { return(threshold.exceedance.duration.index(ci@tmax, ci@annual.factor, ci@jdays, ci@running.pctile.notbase$tx90thresh, ">", spells.can.span.years=spells.can.span.years) * ci@namask.ann$tmax) }

## CSDI
climdex.csdi <- function(ci, spells.can.span.years=FALSE) { return(threshold.exceedance.duration.index(ci@tmin, ci@annual.factor, ci@jdays, ci@running.pctile.notbase$tn10thresh, "<", spells.can.span.years=spells.can.span.years) * ci@namask.ann$tmax) }

## DTR
climdex.dtr <- function(ci, freq=c("monthly", "annual")) { return(mean.daily.temp.range(ci@tmax, ci@tmin, freq.to.factor(ci, match.arg(freq))) * freq.to.namask(ci, match.arg(freq))$tavg) }

## Rx1day
climdex.rx1day <- function(ci, freq=c("monthly", "annual")) { return(nday.consec.prec.max(ci@prec, freq.to.factor(ci, match.arg(freq)), 1) * freq.to.namask(ci, match.arg(freq))$prec) }

## Rx5day
## fclimdex implements Rx5day in a strange way; the running sum series is centered on the last day, and the first day a running sum can be computed for is left out entirely. This results in wet days near a month boundary going into what is arguably the wrong month.
climdex.rx5day <- function(ci, freq=c("monthly", "annual"), center.mean.on.last.day=FALSE) { return(nday.consec.prec.max(ci@prec, freq.to.factor(ci, match.arg(freq)), 5, center.mean.on.last.day) * freq.to.namask(ci, match.arg(freq))$prec) }

## SDII
climdex.sdii <- function(ci) { return(simple.precipitation.intensity.index(ci@prec, ci@annual.factor) * ci@namask.ann$prec) }

## R10mm
climdex.r10mm <- function(ci) { return(number.days.op.threshold(ci@prec, ci@annual.factor, 10, ">=") * ci@namask.ann$prec) }

## R20mm
climdex.r20mm <- function(ci) { return(number.days.op.threshold(ci@prec, ci@annual.factor, 20, ">=") * ci@namask.ann$prec) }

## Rnnmm
climdex.rnnmm <- function(ci, threshold=1) {
  if(!is.numeric(threshold) || length(threshold) != 1) stop("Please specify a single numeric threshold value.");

  return(number.days.op.threshold(ci@prec, ci@annual.factor, threshold, ">=") * ci@namask.ann$prec)
}

## CDD
## Both CDD and CWD in fclimdex do not record the length of consecutive days on transition to a missing value
climdex.cdd <- function(ci, spells.can.span.years=TRUE) { return(spell.length.max(ci@prec, ci@annual.factor, 1, "<", spells.can.span.years) * ci@namask.ann$prec) }

## CWD
climdex.cwd <- function(ci, spells.can.span.years=TRUE) { return(spell.length.max(ci@prec, ci@annual.factor, 1, ">=", spells.can.span.years) * ci@namask.ann$prec) }

## R95pTOT
climdex.r95ptot <- function(ci) { return(total.precip.op.threshold(ci@prec, ci@annual.factor, ci@prec.pctile['r95thresh'], ">") * ci@namask.ann$prec) }

## R99pTOT
climdex.r99ptot <- function(ci) { return(total.precip.op.threshold(ci@prec, ci@annual.factor, ci@prec.pctile['r99thresh'], ">") * ci@namask.ann$prec) }

## PRCPTOT
climdex.prcptot <- function(ci) { return(total.precip.op.threshold(ci@prec, ci@annual.factor, 1, ">=") * ci@namask.ann$prec) }

all.indices <- c('fd', 'su', 'id', 'tr', 'gsl', 'txx', 'tnx', 'txn', 'tnn', 'tn10p', 'tx10p', 'tn90p', 'tx90p', 'wsdi', 'csdi',
                 'dtr', 'rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')

##
## HELPERS FINISHED. IMPLEMENTATIION BELOW.
##

## Input: x: A sequence of booleans
## Input: na.value: The value to replace NA with
## Returns a series of the same length as x; all zeros except at the end of a run of TRUE values, 
## where the value will be the length of the run.
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

## FD, ID, SU, TR, R10mm, R20mm, Rnnmm
## Returns the number of days in each block (as defined by date.factor)
## that exceed (in the case of op=">") or are below (in the case of op=">") the specified threshold.
number.days.op.threshold <- function(temp, date.factor, threshold, op="<") {
  stopifnot(is.numeric(c(temp, threshold)))
  return(tapply.fast(match.fun(op)(temp, threshold), date.factor, sum, na.rm=TRUE))
}

## GSL
## The default mode is meaningless if not annual
## Time series must be contiguous
## NOTE: There is a difference of 1 between our output and fclimdex. See line 637; consider case where start and end day are same. Correct answer is 1 day GSL; their answer is 0 day.
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
    in.gsl <- select.blocks.gt.length(daily.mean.temp >= t.thresh, min.length - 1) | !select.blocks.gt.length(daily.mean.temp < t.thresh, min.length - 1, TRUE)
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

## TN10p, TX10p, TN90p, TX90p
## Computes the number of days for each period (as defined by date.factor) that exceed or are below
## (depending on op) the specified historical threshold for that day.
## Requires use of bootstrap procedure to generate 1961-1990 pctile; see Zhang et al, 2004 
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
  
  return(tapply.fast(dat, date.factor, function(x) { x.nona <- x[!is.na(x)]; if(!length(x.nona)) return(NA); return(sum(x.nona) / length(x.nona) * 100) } ))
}

## WSDI, CSDI
## Computes the total number of days in each period (as defined by date.factor) that exceed or are below
## (depending on op) the specified historical threshold and are part of a spell of at least six days meeting
## the given criteria.
threshold.exceedance.duration.index <- function(daily.temp, date.factor, jdays, thresholds, op=">", min.length=6, spells.can.span.years=TRUE) {
  stopifnot(is.numeric(c(daily.temp, thresholds, min.length)), is.factor(date.factor),
            is.function(match.fun(op)),
            min.length > 0)

  if(spells.can.span.years) {
    periods <- select.blocks.gt.length(match.fun(op)(daily.temp, thresholds[jdays]), min.length - 1)
    return(tapply.fast(periods, date.factor, sum))
  } else {
    ## fclimdex behaviour...
    return(tapply.fast(1:length(daily.temp), date.factor, function(idx) { sum(select.blocks.gt.length(match.fun(op)(daily.temp[idx], thresholds[jdays[idx]]), min.length - 1)) } ))
  }
}

## DTR
## Computes mean diurnal temperature range in each period (as specified by date.factor).
## Max and min temps are assumed to be same length
mean.daily.temp.range <- function(daily.max.temp, daily.min.temp, date.factor) {
  return(tapply.fast(daily.max.temp - daily.min.temp, date.factor, mean))
}

## Rx1day, Rx5day
## Computes the maximum n-day running sum of precipitation in each period (as specified by date.factor).
## Setting center.mean.on.last.day to TRUE reproduces a bug in RX5day in fclimdex, making the results identical
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

## SDII
## Computes the mean precipitation amount in wet days in each period (as specified by date.factor).
## Period for computation of number of wet days shall be the entire range of the data supplied.
simple.precipitation.intensity.index <- function(daily.prec, date.factor) {
  return(tapply.fast(daily.prec, date.factor, function(prec) { idx <- prec >= 1 & !is.na(prec); if(sum(idx) == 0) { return(0); } else { return(sum(prec[idx], na.rm=TRUE) / sum(idx)) } } ))
}

## CDD, CWD
## Computes the longest string of days which exceed or are below (depending on op) the given threshold
## for each period (as specified by date.factor).
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

## R95pTOT, R99pTOT
## Computes the sum of precipitation exceeding or below the specified threshold (depending on op)
## for each period (as specified by date.factor).
total.precip.op.threshold <- function(daily.prec, date.factor, threshold, op) {
  f <- match.fun(op)
  return(tapply.fast(daily.prec, date.factor, function(pr) { return(sum(pr[f(pr, threshold)], na.rm=TRUE)) } ))
}

## Returns an n-day running quantile for each day of data (dimensions c(dpy, q))
## Data is assumed to be padded by floor(n/2) days on either end, and data is assumed to start on the (dpy - floor(n/2) + 1)'th day..
running.quantile <- function(data, n, q, dpy) {
  ret <- .Call("running_quantile_windowed", data, n, q, dpy, DUP=FALSE)
  dim(ret) <- c(length(q), dpy)
  return(t(ret))
}

## Takes an array of booleans; returns an array of booleans where only blocks of TRUE longer than n are still TRUE
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

## Computes quantiles using the C quantile method used by climdex.pcic (equivalent to R's type=8)
climdex.quantile <- function(x, q=c(0, 0.25, 0.5, 0.75, 1)) {
  return(.Call("c_quantile2", as.double(x), q))
}
