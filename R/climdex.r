library(caTools)
library(PCICt)

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

valid.climdexInput <- function(x) {
  errors <- c()

  if(!all(sapply(c("tmax", "tmin", "tavg", "prec", "dates", "annual.factor", "monthly.factor"), function(y) { length(slot(x, y)) }) == length(x@tmax)))
    errors <- c(errors, "Data fields, dates, and factors must all be of the same length")

  ## Check that namask.mon and namask.ann have columns for each of the variables
  if(sum(c("tmax", "tmin", "tavg", "prec") %in% names(x@namask.ann)) != 4)
    errors <- c(errors, "NA mask for annual must contain data for tmax, tmin, tavg, and prec.")

  if(sum(c("tmax", "tmin", "tavg", "prec") %in% names(x@namask.mon)) != 4)
    errors <- c(errors, "NA mask for monthly must contain data for tmax, tmin, tavg, and prec.")

  ## Check that running.pctile.base running.pctile.notbase, and pctile all contain the proper named data.
  bs.pctile.names <- c("tmax10", "tmax90", "tmin10", "tmin90")
  pctile.names <- c("precwet95", "precwet99")

  if(sum(bs.pctile.names %in% names(x@running.pctile.base)) != length(bs.pctile.names))
    errors <- c(errors, "running.pctile.base does not contain at least one of tmax10, tmax90, tmin10, and tmin90.")

  if(sum(bs.pctile.names %in% names(x@running.pctile.notbase)) != length(bs.pctile.names))
    errors <- c(errors, "running.pctile.notbase does not contain at least one of tmax10, tmax90, tmin10, and tmin90.")

  if(sum(pctile.names %in% names(x@pctile)) != length(pctile.names))
    errors <- c(errors, "pctile does not contain at least one of precwet95 and precwet99.")
  
  if(length(errors) == 0)
    return(TRUE)
  else
    return(errors)
}

setClass("climdexInput",
         representation(tmax = "numeric",
                        tmin = "numeric",
                        tavg = "numeric",
                        prec = "numeric",
                        namask.ann = "data.frame",
                        namask.mon = "data.frame",
                        running.pctile.base = "list",
                        running.pctile.notbase = "data.frame",
                        pctile = "numeric",
                        dates = "PCICt",
                        base.range = "PCICt",
                        annual.factor = "factor",
                        monthly.factor = "factor"),
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

create.filled.series <- function(data, data.dates, new.date.sequence) {
  new.data <- rep(NA, length(new.date.sequence))
  data.in.new.data <- (data.dates >= new.date.sequence[1]) & (data.dates <= new.date.sequence[length(new.date.sequence)])
  indices <- round(as.numeric(data.dates[data.in.new.data] - new.date.sequence[1], units="days")) + 1
  new.data[indices] <- data[data.in.new.data]
  return(new.data)
}

get.jdays <- function(dates) {
  return(as.numeric(strftime(dates, "%j", tz="GMT")))
}

get.years <- function(dates) {
  as.numeric(strftime(dates, format="%Y", tz="GMT"))
}

get.jdays.replaced.feb29 <- function(dates) {
  return(unlist(tapply.fast(get.jdays(dates), factor(get.years(dates)), function(x) { if(length(x) == 366) { return(c(1:59, 59, 60:365)) } else { return(x) } })))
}

get.bootstrap.set <- function(dates, bootstrap.range, win.size) {
  bootstrap.win.range <- get.bootstrap.windowed.range(bootstrap.range, win.size)
  return(dates >= bootstrap.win.range[1] & dates <= bootstrap.win.range[2] & strftime(dates, format="%m-%d", tz="GMT") != "02-29")
}

## Input is: bootstrap range, vector of 2 POSIXct, win.size is window size in days (single integer)
get.bootstrap.windowed.range <- function(bootstrap.range, win.size) {
  ## Changed due to a bug in PCICt
  ##window <- as.difftime(floor(win.size / 2), units="days")
  window <- floor(win.size / 2) * 86400
  return(c(bootstrap.range[1] - window, bootstrap.range[2] + window))
}

## Expects POSIXct for all dates
## Do the Zhang boostrapping method described in Xuebin Zhang et al's 2005 paper, "Avoiding Inhomogeneity in Percentile-Based Indices of Temperature Extremes" J.Clim vol 18 pp.1647-1648, "Removing the 'jump'".
zhang.bootstrap.qtile <- function(x, dates, qtiles, bootstrap.range, include.mask=NULL, n=5) {
  window <- floor(n / 2)

  dpy <- ifelse(is.null(attr(dates, "dpy")), 365, attr(dates, "dpy"))
  inset <- get.bootstrap.set(dates, bootstrap.range, n)
  nyears <- floor(sum(inset) / dpy)
  
  if(!is.null(include.mask))
    x[include.mask] <- NA

  bs.data <- x[inset]
  
  ## This routine is written as described in Zhang et al, 2005 as referenced above.
  yidx <- 1:nyears
  d <- sapply(yidx, function(idx.to.omit) {
    ## Index is computed assuming `window` leading days of data...
    omit.index <- window + (1:dpy) + ((idx.to.omit - 1) * dpy)
    return(sapply(yidx[yidx != idx.to.omit], function(idx.to.replace.with) {
      replace.index <- window + (1:dpy) + ((idx.to.replace.with - 1) * dpy)
      bs.data[omit.index] <- bs.data[replace.index]
      return(running.quantile(bs.data, n, qtiles, dpy))
    }))
  } )
  
  dim(d) <- c(dpy, length(qtiles), nyears - 1, nyears)
  d <- aperm(d, perm=c(1, 4, 3, 2))
  ## new dims: dpy, nyears, nyears-1, length(quantiles)

  return(lapply(1:length(qtiles), function(x) { d[,,,x] }))
}

zhang.running.qtile <- function(x, dates, dates.base, qtiles, bootstrap.range, include.mask=NULL, n=5) {
  jdays.idx <- get.jdays.replaced.feb29(dates)
  inset <- get.bootstrap.set(dates.base, bootstrap.range, n)
  dpy <- ifelse(is.null(attr(dates, "dpy")), 365, attr(dates, "dpy"))
  
  if(!is.null(include.mask))
    x[include.mask] <- NA

  bs.data <- x[inset]

  d <- apply(running.quantile(bs.data, n, qtiles, dpy), 2, function(x) { return(x[jdays.idx]) } )

  row.names(d) <- NULL
  return(d)
}

get.na.mask <- function(x, f, threshold) {
  return(c(1, NA)[1 + as.numeric(tapply.fast(is.na(x), f, function(y) { return(sum(y) > threshold) } ))])
}

get.num.days.in.range <- function(x, date.range) {
  return(sum(x >= date.range[1] & x <= date.range[2]))
  
}

climdexInput.raw <- function(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range=c(1961, 1990), n=5) {
  if(length(tmin) != length(tmin.dates))
    stop("Length of tmin data and tmin dates do not match.")
  
  if(length(tmax) != length(tmax.dates))
    stop("Length of tmax data and tmax dates do not match.")
  
  if(length(prec) != length(prec.dates))
    stop("Length of prec data and prec dates do not match.")

  if(!(length(base.range) == 2 && is.numeric(base.range)))
    stop("Invalid base date range; expecting vector of 2 numeric years.")

  if(!(inherits(tmax.dates, "PCICt") && inherits(tmin.dates, "PCICt") && inherits(prec.dates, "PCICt")))
    stop("Dates must all be of class PCICt.")

  if(n != 5)
    warning("Use of n != 5 varies from the Climdex definition. Use at your own risk.")
  
  cal <- attr(tmax.dates, "cal")
  last.day.of.year <- "12-31"
  if(!is.null(attr(tmax.dates, "months")))
    last.day.of.year <- paste("12", attr(tmax.dates, "months")[12], sep="-")
  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)

  ## Check that the base range includes at least a month's worth of data for each variable.
  if(!all(sapply(list(tmax.dates, tmin.dates, prec.dates), get.num.days.in.range, bs.date.range) > 30))
    warning("There is less than a month of data for at least one of tmin, tmax, and prec. Consider revising your base range and/or check your input data.")

  bs.win.date.range <- get.bootstrap.windowed.range(bs.date.range, n)
  all.dates <- c(tmin.dates, tmax.dates, prec.dates, bs.win.date.range)
  date.range <- range(all.dates)
  year.range <- as.numeric(strftime(date.range, "%Y", tz="GMT"))
  new.date.range <- as.PCICt(paste(year.range, c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")
  
  annual.factor <- as.factor(strftime(date.series, "%Y", tz="GMT"))
  monthly.factor <- as.factor(strftime(date.series, "%Y-%m", tz="GMT"))
  
  filled.tmax <- create.filled.series(tmax, tmax.dates, date.series)
  filled.tmin <- create.filled.series(tmin, tmin.dates, date.series)
  filled.prec <- create.filled.series(prec, prec.dates, date.series)
  filled.tavg <- (filled.tmax + filled.tmin) / 2

  filled.list <- list(filled.tmax, filled.tmin, filled.tavg, filled.prec)
  filled.list.names <- c("tmax", "tmin", "tavg", "prec")

  namask.ann <- do.call(data.frame, lapply(filled.list, get.na.mask, annual.factor, 15))
  colnames(namask.ann) <- filled.list.names
  
  namask.mon <- do.call(data.frame, lapply(filled.list, get.na.mask, monthly.factor, 3))
  colnames(namask.mon) <- filled.list.names

  ## DeMorgan's laws FTW
  wet.days <- !(is.na(filled.prec) | filled.prec < 1)

  bs.pctile.base <- do.call(c, lapply(filled.list[1:2], zhang.bootstrap.qtile, date.series, c(0.1, 0.9), bs.date.range, n=n))
  bs.pctile <- do.call(data.frame, lapply(filled.list[1:2], zhang.running.qtile, date.series, date.series, c(0.1, 0.9), bs.date.range, n=n))

  inset <- date.series >= new.date.range[1] & date.series <= new.date.range[2] & !is.na(filled.prec) & wet.days
  pctile <- quantile(filled.prec[inset], c(0.95, 0.99))
  
  names(bs.pctile.base) <- c("tmax10", "tmax90", "tmin10", "tmin90")
  names(bs.pctile) <- c("tmax10", "tmax90", "tmin10", "tmin90")
  names(pctile) <- c("precwet95", "precwet99")
  
  return(new("climdexInput", tmax=filled.tmax, tmin=filled.tmin, tavg=filled.tavg, prec=filled.prec, namask.ann=namask.ann, namask.mon=namask.mon, running.pctile.base=bs.pctile.base, running.pctile.notbase=bs.pctile, pctile=pctile, dates=date.series, base.range=bs.date.range, annual.factor=annual.factor, monthly.factor=monthly.factor))
}

climdexInput.csv <- function(tmax.file, tmin.file, prec.file, data.columns=list(tmin="tmin", tmax="tmax", prec="prec"), base.range=c(1961, 1990), na.strings=NULL, cal="gregorian", date.types=NULL, n=5) {
  if(sum(c("tmax", "tmin", "prec") %in% names(data.columns)) != 3)
    stop("Must specify names of all data columns (tmin, tmax, prec).")

  if(missing(date.types))
    date.types <- list(list(fields=c("year", "jday"), format="%Y %j"),
                       list(fields=c("year", "month", "day"), format="%Y %m %d"))
  else
    if(any(!sapply(date.types, function(x) { return(sum(c("fields", "format") %in% names(x)) == 2 && is.character(x$fields) && is.character(x$format)) } )))
      stop("Invalid date.types specified. See ?climdexInput.csv .")
  
  tmin.dat <- read.csv(tmin.file, na.strings=na.strings)
  tmax.dat <- read.csv(tmax.file, na.strings=na.strings)
  prec.dat <- read.csv(prec.file, na.strings=na.strings)

  ## This is to deal with fclimdex's broken input data, which includes February 31st, amongst other joyous things
  tmin.dat <- tmin.dat[!is.na(tmin.dat[,4]),]
  tmax.dat <- tmax.dat[!is.na(tmax.dat[,4]),]
  prec.dat <- prec.dat[!is.na(prec.dat[,4]),]
  
  if(!(data.columns$tmin %in% names(tmin.dat) & data.columns$tmax %in% names(tmax.dat) & data.columns$prec %in% names(prec.dat)))
    stop("Data columns not found in data.")
  
  tmin.dates <- get.date.field(tmin.dat, cal, date.types)
  tmax.dates <- get.date.field(tmax.dat, cal, date.types)
  prec.dates <- get.date.field(prec.dat, cal, date.types)

  return(climdexInput.raw(tmax.dat[,data.columns$tmax], tmin.dat[,data.columns$tmin], prec.dat[,data.columns$prec], tmax.dates, tmin.dates, prec.dates, base.range))
}

## Temperature units: degrees C
## Precipitation units: mm per unit time

## Status:
## FD: Annual. Differences of 1-4 days in some years.
climdex.fd <- function(ci) { return(number.days.op.threshold(ci@tmin, ci@annual.factor, 0, "<") * ci@namask.ann$tmin) }

## SU: Annual. Differences of 1-4 days in some years.
climdex.su <- function(ci) { return(number.days.op.threshold(ci@tmax, ci@annual.factor, 25, ">") * ci@namask.ann$tmax) }

## ID: Annual. Differences of 1-4 days in some years.
climdex.id <- function(ci) { return(number.days.op.threshold(ci@tmax, ci@annual.factor, 0, "<") * ci@namask.ann$tmax) }

## TR: Annual. Differences of 1-4 days in some years.
climdex.tr <- function(ci) { return(number.days.op.threshold(ci@tmin, ci@annual.factor, 20, ">") * ci@namask.ann$tmin) }

## GSL: Annual. Should work, needs more testing; is imprecise around date of Jul 1. Creates GSL 1 day longer than fclimdex due to off-by-one in fclimdex.
climdex.gsl <- function(ci) { return(growing.season.length(ci@tavg, ci@annual.factor) * ci@namask.ann$tavg) }

## TXx: Monthly. Exact match.
climdex.txx <- function(ci) { return(tapply.fast(ci@tmax, ci@monthly.factor, max) * ci@namask.mon$tmax) }

## TNx: Monthly. Exact match.
climdex.tnx <- function(ci) { return(tapply.fast(ci@tmin, ci@monthly.factor, max) * ci@namask.mon$tmin) }

## TXn: Monthly. Exact match.
climdex.txn <- function(ci) { return(tapply.fast(ci@tmax, ci@monthly.factor, min) * ci@namask.mon$tmax) }

## TNn: Monthly. Exact match.
climdex.tnn <- function(ci) { return(tapply.fast(ci@tmin, ci@monthly.factor, min) * ci@namask.mon$tmin) }

## TN10p: Monthly. Pattern matches, but still significant differences.
## Our implementation currently follows the example set by fclimdex for dealing with missing values, which is wrong; it biases results upwards when missing values are present.
#percent.days.op.threshold <- function(temp, dates, annual.factor, threshold.outside.base, base.thresholds, base.range, op='<') {

climdex.tn10p <- function(ci) { return(percent.days.op.threshold(ci@tmin, ci@dates, ci@monthly.factor, ci@running.pctile.notbase$tmin10, ci@running.pctile.base$tmin10, ci@base.range, "<") * ci@namask.mon$tmin) }

## TX10p: Monthly. Pattern matches, but still significant differences.
climdex.tx10p <- function(ci) { return(percent.days.op.threshold(ci@tmax, ci@dates, ci@monthly.factor, ci@running.pctile.notbase$tmax10, ci@running.pctile.base$tmax10, ci@base.range, "<") * ci@namask.mon$tmax) }

## TN90p: Monthly. Pattern matches, but still significant differences.
climdex.tn90p <- function(ci) { return(percent.days.op.threshold(ci@tmin, ci@dates, ci@monthly.factor, ci@running.pctile.notbase$tmin90, ci@running.pctile.base$tmin90, ci@base.range, ">") * ci@namask.mon$tmin) }

## TX90p: Monthly. Pattern matches, but still significant differences.
climdex.tx90p <- function(ci) { return(percent.days.op.threshold(ci@tmax, ci@dates, ci@monthly.factor, ci@running.pctile.notbase$tmax90, ci@running.pctile.base$tmax90, ci@base.range, ">") * ci@namask.mon$tmax) }

## WSDI: Annual. Significant differences.
## fclimdex implements CSDI and WSDI incorrectly; it adds the entire spell to the year in which the spell ended. This code sums up the days which were part of the spell.
climdex.wsdi <- function(ci) { return(threshold.exceedance.duration.index(ci@tmax, ci@annual.factor, ci@running.pctile.notbase$tmax90, ">") * ci@namask.ann$tavg) }

## CSDI: Annual. Pattern matches but significant differences exist.
climdex.csdi <- function(ci) { return(threshold.exceedance.duration.index(ci@tmin, ci@annual.factor, ci@running.pctile.notbase$tmin10, "<") * ci@namask.ann$tavg) }

## DTR: Monthly. Differences in some samples at the 3rd decimal place.
climdex.dtr <- function(ci) { return(mean.daily.temp.range(ci@tmax, ci@tmin, ci@monthly.factor) * ci@namask.mon$tavg) }

## Rx1day: Monthly. Exact match.
climdex.rx1day <- function(ci) { return(max.nday.consec.prec(ci@prec, ci@monthly.factor, 1) * ci@namask.mon$prec) }

## Rx5day: Monthly. Code should be correct.
## fclimdex implements Rx5day incorrectly; the running sum series is off by 2 days, and the first day a running sum can be computed for is left out entirely. This results in wet days near a month boundary going into a different month.
climdex.rx5day <- function(ci) { return(max.nday.consec.prec(ci@prec, ci@monthly.factor, 5) * ci@namask.mon$prec) }

## SDII: Annual. Small differences due to fclimdex's rounding to 1 decimal place.
climdex.sdii <- function(ci) { return(simple.precipitation.intensity.index(ci@prec, ci@annual.factor) * ci@namask.ann$prec) }

## R10mm: Annual. Exact match.
climdex.r10mm <- function(ci) { return(number.days.op.threshold(ci@prec, ci@annual.factor, 10, ">=") * ci@namask.ann$prec) }

## R20mm: Annual. Exact match.
climdex.r20mm <- function(ci) { return(number.days.op.threshold(ci@prec, ci@annual.factor, 20, ">=") * ci@namask.ann$prec) }

## Rnnmm: Annual. Exact match.
climdex.rnnmm <- function(ci, threshold) {
  if(!is.numeric(threshold) || length(threshold) != 1) stop("Please specify a single numeric threshold value.");

  return(number.days.op.threshold(ci@prec, ci@annual.factor, threshold, ">=") * ci@namask.ann$prec)
}

## Both CDD and CWD in fclimdex do not record the length of consecutive days on transition to a missing value
## CDD: Annual. Exact match.
climdex.cdd <- function(ci) { return(max.length.spell(ci@prec, ci@annual.factor, 1, "<") * ci@namask.ann$prec) }

## CWD: Annual. Exact match.
climdex.cwd <- function(ci) { return(max.length.spell(ci@prec, ci@annual.factor, 1, ">=") * ci@namask.ann$prec) }

## R95pTOT: Annual. Exact match.
climdex.r95ptot <- function(ci) { return(total.precip.op.threshold(ci@prec, ci@annual.factor, ci@pctile['precwet95'], ">") * ci@namask.ann$prec) }

## R99pTOT: Annual. Exact match.
climdex.r99ptot <- function(ci) { return(total.precip.op.threshold(ci@prec, ci@annual.factor, ci@pctile['precwet99'], ">") * ci@namask.ann$prec) }

## PRCPTOT: Annual. Exact match.
climdex.prcptot <- function(ci) { return(total.precip.op.threshold(ci@prec, ci@annual.factor, 1, ">=") * ci@namask.ann$prec) }

all.indicies <- c('fd', 'su', 'id', 'tr', 'gsl', 'txx', 'tnx', 'txn', 'tnn', 'tn10p', 'tx10p', 'tn90p', 'tx90p', 'wsdi', 'csdi',
                  'dtr', 'rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')

##
## HELPERS FINISHED. IMPLEMENTATIION BELOW.
##


get.series.lengths.at.ends <- function(x) {
  n <- length(x)
  if(n == 1)
    return(as.numeric(x))

  res <- rep(0, n)
  x[is.na(x)] <- FALSE

  ## Compare series to lag-1 and lag+1 series; false added to trigger state transition from TRUE at ends of series
  start <- which(x & !(c(FALSE, x[1:(n - 1)])))
  end <- which(x & !(c(x[2:n], FALSE)))
  res[end] <- end - start + 1
  return(res)
}

## FD, ID, SU, TR, R10mm, R20mm, Rnnmm
number.days.op.threshold <- function(temp, date.factor, threshold, op="<") {
  stopifnot(is.numeric(c(temp, threshold)))
  return(tapply.fast(match.fun(op)(temp, threshold), date.factor, sum, na.rm=TRUE))
}

## GSL
## Meaningless if not annual
## Time series must be contiguous
## NOTE: There is a difference of 1 between our output and fclimdex. See line 637; consider case where start and end day are same. Correct answer is 1 day GSL; their answer is 0 day.
growing.season.length <- function(daily.mean.temp, date.factor,
                                  min.length=6, t.thresh=5) {
  return(tapply.fast(daily.mean.temp, date.factor, function(ts) {
    ts.len<- length(ts)
    ts.mid <- floor(ts.len / 2)
    gs.begin <- which(select.blocks.gt.length(ts[1:(ts.mid-1)] > t.thresh, min.length - 1))

    ## Growing season actually ends the day -before- the sequence of sketchy days
    gs.end <- which(select.blocks.gt.length(ts[ts.mid:ts.len] < t.thresh, min.length - 1)) - 1

    if(length(gs.begin) == 0) {
      return(0)
    } else if(length(gs.end) == 0) {
      return(ts.len - gs.begin[1] + 1)
    } else {
      return(gs.end[1] - gs.begin[1] + 1 + ts.mid)
    }
  } ))
}

## TN10p, TX10p, TN90p, TX90p
## Requires use of bootstrap procedure to generate 1961-1990 pctile; see Zhang et al, 2004 
percent.days.op.threshold <- function(temp, dates, date.factor, threshold.outside.base, base.thresholds, base.range, op='<') {
  f <- match.fun(op)
  
  dat <- f(temp, threshold.outside.base)
  inset <- dates >= base.range[1] & dates <= base.range[2]
  
  if(sum(inset) > 0) {
    years.base <- get.years(dates[inset])
    jdays.base <- get.jdays.replaced.feb29(dates[inset])
    
    temp.base <- temp[inset]
    years.base.range <- range(years.base)
    byrs <- (years.base.range[2] - years.base.range[1] + 1)
    year.base.list <- years.base.range[1]:years.base.range[2]

    d <- lapply(1:byrs, function(x) { yset <-  which(years.base == year.base.list[x]); sapply(1:(byrs - 1), function(y) { f(temp.base[yset], base.thresholds[jdays.base[yset],x,y]) } ) })
    ## This repeats a bug (or at least, debatable decision) in fclimdex where they always divide by byrs - 1 even when there are NAs (missing values) in the thresholds
    dat[inset] <- unlist(lapply(d, apply, 1, function(x) { sum(as.numeric(x), na.rm=TRUE) } ) ) / (byrs - 1)
  }
  
  return(tapply.fast(dat, date.factor, function(x) { x.nona <- x[!is.na(x)]; if(!length(x.nona)) return(NA); return(sum(x.nona) / length(x.nona) * 100) } ))
}

## WSDI, CSDI
## Thresholds appear to be for each block of 5 days of a year...
threshold.exceedance.duration.index <- function(daily.max.temp, date.factor, warm.thresholds, op=">", min.length=6) {
  stopifnot(is.numeric(c(daily.max.temp, warm.thresholds, min.length)), is.factor(date.factor),
            is.function(match.fun(op)),
            min.length > 0)

  periods <- select.blocks.gt.length(match.fun(op)(daily.max.temp, warm.thresholds), min.length - 1)
  return(tapply.fast(periods, date.factor, sum))
}

## DTR
## Max and min temps are assumed to be same length
mean.daily.temp.range <- function(daily.max.temp, daily.min.temp, date.factor) {
  return(tapply.fast(daily.max.temp - daily.min.temp, date.factor, mean))
}

## Rx1day, Rx5day
max.nday.consec.prec <- function(daily.prec, date.factor, ndays) {
  if(ndays == 1) {
    return(tapply.fast(daily.prec, date.factor, max))
  } else {
    ## Ends of the data will be de-emphasized (padded with zero precip data); NAs replaced with 0
    new.series <- c(rep(0, floor(ndays / 2)), daily.prec, rep(0, floor(ndays / 2)))
    new.series[is.na(new.series)] <- 0
    prec.runsum <- runmean(new.series, k=ndays, endrule="trim") * ndays
    ## Uncommenting this introduces the bug in RX5day in fclimdex, making the results identical
    ##prec.runsum <- c(0, 0, prec.runsum[1:(length(prec.runsum) - floor(ndays / 2))])
    return(tapply.fast(prec.runsum, date.factor, max))
  }
}

## SDII
## Period for computation of number of wet days shall be the entire range of the data supplied.
simple.precipitation.intensity.index <- function(daily.prec, date.factor) {
  return(tapply.fast(daily.prec, date.factor, function(prec) { idx <- prec >= 1 & !is.na(prec); if(sum(idx) == 0) { return(0); } else { return(sum(prec[idx], na.rm=TRUE) / sum(idx)) } } ))
}

## CDD, CWD
max.length.spell <- function(daily.prec, date.factor, threshold, op) {
  bools <- match.fun(op)(daily.prec, threshold)
  return(tapply.fast(get.series.lengths.at.ends(bools), date.factor, function(x) { return(max(x)) } ))
}

## R95pTOT, R99pTOT
total.precip.op.threshold <- function(daily.prec, date.factor, threshold, op) {
  f <- match.fun(op)
  return(tapply.fast(1:length(daily.prec), date.factor, function(x) { return(sum(daily.prec[x[f(daily.prec[x], threshold)]], na.rm=TRUE)) } ))
}

## Returns an n-day running quantile for each day of data
## Data is assumed to be padded by floor(n/2) days on either end, and data is assumed to start on the (dpy - floor(n/2) + 1)'th day..
running.quantile <- function(data, n, q, dpy) {
  ret <- .Call("running_quantile_windowed", data, n, q, dpy, DUP=FALSE)
  dim(ret) <- c(length(q), dpy)
  ##browser()
  return(t(ret))
}

## Returns an n-day running quantile for each day of data
## Data is assumed to be 365 days per year, data is assumed to be padded by floor(n/2) days on either end, and data is assumed to start on the (365 - floor(n/2) + 1)'th day..
running.quantile.r <- function(data, n, q, include.mask=NULL) {
  window <- floor(n / 2)
  true.data.length <- length(data) - 2 * window

  ## Doesn't seem to make a difference (this is behaviour of fclimdex)
  ##data[1:window] <- data[window + 1]
  ##data[length(data) - (0:(window - 1))] <- data[length(data) - window]
  
  ## Apply include mask
  if(!is.null(include.mask))
    data[include.mask] <- NA

  ## Create n lists of indices 
  data.perm <- unlist(lapply(1:n, function(x) { return(data[x:(true.data.length + x - 1)]) }))
  dim(data.perm) <- c(365, ceiling(length(data.perm) / 365))

  ## Transpose so apply goes faster
  data.perm <- t(data.perm)

  ## Return transposed data so that major dim is proportional to length of q
  d <- t(apply(data.perm, 2, quantile, q, na.rm=TRUE, type=8))
  return(d)
}

## Takes an array of booleans; returns an array of booleans where only blocks of TRUE longer than n are still TRUE
select.blocks.gt.length <- function(d, n) {
  stopifnot(is.logical(d), is.numeric(n))

  if(n <= 1)
    return(d)

  if(n >= length(d))
    return(rep(FALSE, length(d)))

  d[is.na(d)] <- FALSE
  
  d2 <- Reduce(function(x, y) { return(c(rep(FALSE, y), d[1:(length(d) - y)]) & x) }, 1:n, d)
  return(Reduce(function(x, y) { return(c(d2[(y + 1):length(d2)], rep(FALSE, y)) | x) }, 1:n, d2))
}
