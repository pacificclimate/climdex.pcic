library(caTools)
library(PCICt)

temp.quantiles <- c(10, 90)
prec.quantiles <- c(95, 99)

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
    d <- .Call("running_quantile_windowed_bootstrap", bs.data, n, qtiles, dpy, DUP=FALSE)
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

## Get quantiles for out-of-base period; for use when computing indices on future data using historical quantiles.
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

## Creates climdexInput data structure given input.
## Temp. units: degrees C, Prec. units: mm per day
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

## Creates climdexInput data structure given input CSV files consisting of date columns and data.
## Temp. units: degrees C, Prec. units: mm per day
climdexInput.csv <- function(tmax.file=NULL, tmin.file=NULL, prec.file=NULL,
                             data.columns=list(tmin="tmin", tmax="tmax", prec="prec"), base.range=c(1961, 1990),
                             na.strings=NULL, cal="gregorian", date.types=NULL, n=5, northern.hemisphere=TRUE,
                             pad.data.with.first.last.values=FALSE, tavg.file=NULL, quantiles=NULL) {
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
  
  return(climdexInput.raw(tmax=tmax$dat, tmin=tmin$dat, prec=prec$dat, tmax.dates=tmax$dates, tmin.dates=tmin$dates, prec.dates=prec$dates, base.range=base.range, n=n, northern.hemisphere=northern.hemisphere, pad.data.with.first.last.values=pad.data.with.first.last.values, tavg=tavg$dat, tavg.dates=tavg$dates, quantiles=quantiles))
}
  
## Returns appropriate factor for given time frequency
freq.to.factor <- function(ci, freq) { switch(freq, annual=ci@annual.factor, monthly=ci@monthly.factor, NULL) }

## Returns appropriate NA mask for given time frequency
freq.to.namask <- function(ci, freq) { switch(freq, annual=ci@namask.ann, monthly=ci@namask.mon, NULL) }

## FD: Differences of 1-4 days in some years.
climdex.fd <- function(ci) { stopifnot(!is.null(ci@data$tmin)); return(number.days.op.threshold(ci@data$tmin, ci@annual.factor, 0, "<") * ci@namask.ann$tmin) }

## SU: Differences of 1-4 days in some years.
climdex.su <- function(ci) { stopifnot(!is.null(ci@data$tmax)); return(number.days.op.threshold(ci@data$tmax, ci@annual.factor, 25, ">") * ci@namask.ann$tmax) }

## ID: Differences of 1-4 days in some years.
climdex.id <- function(ci) { stopifnot(!is.null(ci@data$tmax)); return(number.days.op.threshold(ci@data$tmax, ci@annual.factor, 0, "<") * ci@namask.ann$tmax) }

## TR: Differences of 1-4 days in some years.
climdex.tr <- function(ci) { stopifnot(!is.null(ci@data$tmin)); return(number.days.op.threshold(ci@data$tmin, ci@annual.factor, 20, ">") * ci@namask.ann$tmin) }

## GSL: Modes other than "GSL" may not be well tested.
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

## TXx
climdex.txx <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax)); return(tapply.fast(ci@data$tmax, freq.to.factor(ci, match.arg(freq)), max) * freq.to.namask(ci, match.arg(freq))$tmax) }

## TNx
climdex.tnx <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin)); return(tapply.fast(ci@data$tmin, freq.to.factor(ci, match.arg(freq)), max) * freq.to.namask(ci, match.arg(freq))$tmin) }

## TXn
climdex.txn <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax)); return(tapply.fast(ci@data$tmax, freq.to.factor(ci, match.arg(freq)), min) * freq.to.namask(ci, match.arg(freq))$tmax) }

## TNn
climdex.tnn <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin)); return(tapply.fast(ci@data$tmin, freq.to.factor(ci, match.arg(freq)), min) * freq.to.namask(ci, match.arg(freq))$tmin) }

## TN10p
## Our implementation currently follows the example set by fclimdex for dealing with missing values, which is wrong; it biases results upwards when missing values are present.
climdex.tn10p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(percent.days.op.threshold(ci@data$tmin, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmin$outbase$q10, ci@quantiles$tmin$inbase$q10, ci@base.range, "<") * freq.to.namask(ci, match.arg(freq))$tmin) }

## TX10p
climdex.tx10p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmax$outbase$q10, ci@quantiles$tmax$inbase$q10, ci@base.range, "<") * freq.to.namask(ci, match.arg(freq))$tmax) }

## TN90p
climdex.tn90p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(percent.days.op.threshold(ci@data$tmin, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmin$outbase$q90, ci@quantiles$tmin$inbase$q90, ci@base.range, ">") * freq.to.namask(ci, match.arg(freq))$tmin) }

## TX90p
climdex.tx90p <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, freq.to.factor(ci, match.arg(freq)), ci@quantiles$tmax$outbase$q90, ci@quantiles$tmax$inbase$q90, ci@base.range, ">") * freq.to.namask(ci, match.arg(freq))$tmax) }

## WSDI
climdex.wsdi <- function(ci, spells.can.span.years=FALSE) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(threshold.exceedance.duration.index(ci@data$tmax, ci@annual.factor, ci@jdays, ci@quantiles$tmax$outbase$q90, ">", spells.can.span.years=spells.can.span.years) * ci@namask.ann$tmax) }

## CSDI
climdex.csdi <- function(ci, spells.can.span.years=FALSE) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(threshold.exceedance.duration.index(ci@data$tmin, ci@annual.factor, ci@jdays, ci@quantiles$tmin$outbase$q10, "<", spells.can.span.years=spells.can.span.years) * ci@namask.ann$tmax) }

## DTR
climdex.dtr <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@data$tmax) && !is.null(ci@data$tavg)); return(mean.daily.temp.range(ci@data$tmax, ci@data$tmin, freq.to.factor(ci, match.arg(freq))) * freq.to.namask(ci, match.arg(freq))$tavg) }

## Rx1day
climdex.rx1day <- function(ci, freq=c("monthly", "annual")) { stopifnot(!is.null(ci@data$prec)); return(nday.consec.prec.max(ci@data$prec, freq.to.factor(ci, match.arg(freq)), 1) * freq.to.namask(ci, match.arg(freq))$prec) }

## Rx5day
## fclimdex implements Rx5day in a strange way; the running sum series is centered on the last day, and the first day a running sum can be computed for is left out entirely. This results in wet days near a month boundary going into what is arguably the wrong month.
climdex.rx5day <- function(ci, freq=c("monthly", "annual"), center.mean.on.last.day=FALSE) { stopifnot(!is.null(ci@data$prec)); return(nday.consec.prec.max(ci@data$prec, freq.to.factor(ci, match.arg(freq)), 5, center.mean.on.last.day) * freq.to.namask(ci, match.arg(freq))$prec) }

## SDII
climdex.sdii <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(simple.precipitation.intensity.index(ci@data$prec, ci@annual.factor) * ci@namask.ann$prec) }

## R10mm
climdex.r10mm <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(number.days.op.threshold(ci@data$prec, ci@annual.factor, 10, ">=") * ci@namask.ann$prec) }

## R20mm
climdex.r20mm <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(number.days.op.threshold(ci@data$prec, ci@annual.factor, 20, ">=") * ci@namask.ann$prec) }

## Rnnmm
climdex.rnnmm <- function(ci, threshold=1) {
  stopifnot(!is.null(ci@data$prec));
  if(!is.numeric(threshold) || length(threshold) != 1) stop("Please specify a single numeric threshold value.");

  return(number.days.op.threshold(ci@data$prec, ci@annual.factor, threshold, ">=") * ci@namask.ann$prec)
}

## CDD
## Both CDD and CWD in fclimdex do not record the length of consecutive days on transition to a missing value
climdex.cdd <- function(ci, spells.can.span.years=TRUE) { stopifnot(!is.null(ci@data$prec)); return(spell.length.max(ci@data$prec, ci@annual.factor, 1, "<", spells.can.span.years) * ci@namask.ann$prec) }

## CWD
climdex.cwd <- function(ci, spells.can.span.years=TRUE) { stopifnot(!is.null(ci@data$prec)); return(spell.length.max(ci@data$prec, ci@annual.factor, 1, ">=", spells.can.span.years) * ci@namask.ann$prec) }

## R95pTOT
climdex.r95ptot <- function(ci) { stopifnot(!is.null(ci@data$prec) && !is.null(ci@quantiles$prec)); return(total.precip.op.threshold(ci@data$prec, ci@annual.factor, ci@quantiles$prec['q95'], ">") * ci@namask.ann$prec) }

## R99pTOT
climdex.r99ptot <- function(ci) { stopifnot(!is.null(ci@data$prec) && !is.null(ci@quantiles$prec)); return(total.precip.op.threshold(ci@data$prec, ci@annual.factor, ci@quantiles$prec['q99'], ">") * ci@namask.ann$prec) }

## PRCPTOT
climdex.prcptot <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(total.precip.op.threshold(ci@data$prec, ci@annual.factor, 1, ">=") * ci@namask.ann$prec) }

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
  
  ret <- tapply.fast(dat, date.factor, function(x) { return(sum(x, na.rm=TRUE) / sum(!is.na(x))); }) * 100
  ret[is.nan(ret)] <- NA
  return(ret)
}

## WSDI, CSDI
## Computes the total number of days in each period (as defined by date.factor) that exceed or are below
## (depending on op) the specified historical threshold and are part of a spell of at least six days meeting
## the given criteria.
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
