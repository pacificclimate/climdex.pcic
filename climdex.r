library(caTools)
library(abind)

setClass("climdexInput",
         representation(tmax = "numeric",
                        tmin = "numeric",
                        tavg = "numeric",
                        prec = "numeric",
                        namask.ann = "data.frame",
                        namask.mon = "data.frame",
                        bs.pctile = "data.frame",
                        annual.factor = "factor",
                        monthly.factor = "factor")
         )


## Returns POSIXct field or dies
get.date.field <- function(input.data) {
  date.types <- list(list(fields=c("year", "jday"), format="%Y %j"),
                     list(fields=c("year", "month", "day"), format="%Y %m %d"))
  valid.date.types <- sapply(date.types, function(x) { return(!inherits(try(input.data[,x$fields], silent=TRUE), "try-error")) })

  if(sum(valid.date.types) == 0) {
    stop("Could not find a workable set of date fields")
  }

  date.type <- date.types[[which(valid.date.types)[1]]]
  date.strings <- do.call(paste, input.data[,date.type$fields])
  return(as.POSIXct(date.strings, format=date.type$format, tz="GMT"))
}

create.filled.series <- function(data, data.dates, new.date.sequence) {
  new.data <- rep(NA, length(new.date.sequence))
  data.in.new.data <- data.dates >= new.date.sequence[1] & data.dates <= new.date.sequence[length(new.date.sequence)]
  indices <- round(as.numeric(data.dates[data.in.new.data] - new.date.sequence[1], units="days")) + 1
  new.data[indices] <- data[data.in.new.data]
  return(new.data)
}

## Expects POSIXct for all dates
## Do the Zhang boostrapping method described in Xuebin Zhang et al's 2005 paper, "Avoiding Inhomogeneity in Percentile-Based Indices of Temperature Extremes" J.Clim vol 18 pp.1647-1648, "Removing the 'jump'".
## Except don't, because that's not what the Fortran code does.
zhang.bootstrap.qtile <- function(x, dates, qtiles, bootstrap.range, include.mask=NULL) {
  jdays.all <- as.numeric(strftime(dates, "%j", tz="GMT"))
  inset <- dates >= bootstrap.range[1] & dates <= bootstrap.range[2] & !is.na(x)
  years <- as.numeric(strftime(dates[inset], format="%Y", tz="GMT"))
  bs.data <- x[inset]
  jdays <- jdays.all[inset]
  year.list <- unique(years)

  if(!is.null(include.mask))
    include.mask <- include.mask[inset]
  
  ## This routine is written as described in Zhang et al, 2005 as referenced above. However, the Fortran code simply doesn't use this method.
  ##omit.year.data <- do.call(abind, c(lapply(year.list[1:(length(year.list) - 1)], function(year.to.omit) {
  ##  year.to.repeat <- year.to.omit + 1
  ##  my.set <- c(which(!(years == year.to.omit)), which(years == year.to.repeat))
  ##  return(running.quantile(bs.data[my.set], jdays[my.set], 5, qtiles))
  ##} ), along=3))
  ##return(apply(omit.year.data, c(1, 2), mean))

  d <- apply(running.quantile(bs.data, jdays, 5, qtiles, include.mask), 2, function(x) { return(x[jdays.all]) } )
  row.names(d) <- NULL
  return(d)
}

get.na.mask <- function(x, f, threshold) {
  return(c(1, NA)[1 + as.numeric(tapply(is.na(x), f, function(y) { return(sum(y) > threshold) } ))])
}

climdexInput <- function(tmax.file, tmin.file, prec.file, data.columns=list(tmin="tmin", tmax="tmax", prec="prec"), base.range=c(1961, 1990), pctile=c(10, 90)) {
  tmin.dat <- read.csv(tmin.file)
  tmax.dat <- read.csv(tmax.file)
  prec.dat <- read.csv(prec.file)

  if(!(data.columns$tmin %in% names(tmin.dat) & data.columns$tmax %in% names(tmax.dat) & data.columns$prec %in% names(prec.dat))) {
    stop("Data columns not found in data.")
  }
  
  tmin.dates <- get.date.field(tmin.dat)
  tmax.dates <- get.date.field(tmax.dat)
  prec.dates <- get.date.field(prec.dat)

  date.range <- range(c(tmin.dates, tmax.dates, prec.dates))
  year.range <- as.numeric(strftime(date.range, "%Y"))
  new.date.range <- as.POSIXct(paste(year.range, c("01-01", "12-31"), sep="-"), tz="GMT")
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")

  annual.factor <- as.factor(strftime(date.series, "%Y", tz="GMT"))
  monthly.factor <- as.factor(strftime(date.series, "%Y-%m", tz="GMT"))

  filled.tmax <- create.filled.series(tmax.dat[,data.columns$tmax], tmax.dates, date.series)
  filled.tmin <- create.filled.series(tmin.dat[,data.columns$tmin], tmin.dates, date.series)
  filled.tavg <- (filled.tmax + filled.tmin) / 2
  filled.prec <- create.filled.series(prec.dat[,data.columns$prec], prec.dates, date.series)

  filled.list <- list(filled.tmax, filled.tmin, filled.tavg, filled.prec)
  filled.list.names <- c("tmax", "tmin", "tavg", "prec")
  
  namask.ann <- do.call(data.frame, lapply(filled.list, get.na.mask, annual.factor, 15))
  colnames(namask.ann) <- filled.list.names
  
  namask.mon <- do.call(data.frame, lapply(filled.list, get.na.mask, monthly.factor, 3))
  colnames(namask.mon) <- filled.list.names

  ## DeMorgan's laws FTW
  wet.days <- !(is.na(filled.prec) | filled.prec < 1)

  bs.pctile <- do.call(data.frame, c(lapply(filled.list[1:2], zhang.bootstrap.qtile, date.series, c(0.1, 0.9), new.date.range), list(zhang.bootstrap.qtile(filled.prec, date.series, c(0.95, 0.99), new.date.range, wet.days))))
  ##browser()
  
  names(bs.pctile) <- c("tmax10", "tmax90", "tmin10", "tmin90", "precwet95", "precwet99")

  return(new("climdexInput", tmax=filled.tmax, tmin=filled.tmin, tavg=filled.tavg, prec=filled.prec, namask.ann=namask.ann, namask.mon=namask.mon, bs.pctile=bs.pctile, annual.factor=annual.factor, monthly.factor=monthly.factor))
}

## Temperature units: degrees C
## Precipitation units: mm per unit time

## Status:
## FD: Annual. Done
climdex.fd <- function(ci) { return(number.days.op.threshold(ci@tmin, ci@annual.factor, 0, "<") * ci@namask.ann$tmin) }

## SU: Annual. Done
climdex.su <- function(ci) { return(number.days.op.threshold(ci@tmax, ci@annual.factor, 25, ">") * ci@namask.ann$tmax) }

## ID: Annual. Done
climdex.id <- function(ci) { return(number.days.op.threshold(ci@tmax, ci@annual.factor, 0, "<") * ci@namask.ann$tmax) }

## TR: Annual. Done
climdex.tr <- function(ci) { return(number.days.op.threshold(ci@tmin, ci@annual.factor, 20, ">") * ci@namask.ann$tmin) }

## GSL: Annual. Should work, needs more testing; is imprecise around date of Jul 1
climdex.gsl <- function(ci) { return(growing.season.length(ci@tavg, ci@annual.factor) * ci@namask.ann$tavg) }

## TXx: Monthly. Done
climdex.txx <- function(ci) { return(tapply(ci@tmax, ci@monthly.factor, max) * ci@namask.mon$tmax) }

## TNx: Monthly. Done
climdex.tnx <- function(ci) { return(tapply(ci@tmin, ci@monthly.factor, max) * ci@namask.mon$tmin) }

## TXn: Monthly. Done
climdex.txn <- function(ci) { return(tapply(ci@tmax, ci@monthly.factor, min) * ci@namask.mon$tmax) }

## TNn: Monthly. Done
climdex.tnn <- function(ci) { return(tapply(ci@tmin, ci@monthly.factor, min) * ci@namask.mon$tmin) }

## TN10p: Monthly.
climdex.tn10p <- function(ci) { return(percent.days.op.threshold(ci@tmin, ci@monthly.factor, ci@bs.pctile$tmin10, "<") * ci@namask.mon$tmin) }

## TX10p: Monthly.
climdex.tx10p <- function(ci) { return(percent.days.op.threshold(ci@tmax, ci@monthly.factor, ci@bs.pctile$tmax10, "<") * ci@namask.mon$tmax) }

## TN90p: Monthly.
climdex.tn90p <- function(ci) { return(percent.days.op.threshold(ci@tmin, ci@monthly.factor, ci@bs.pctile$tmin90, ">") * ci@namask.mon$tmin) }

## TX90p: Monthly.
climdex.tx90p <- function(ci) { return(percent.days.op.threshold(ci@tmax, ci@monthly.factor, ci@bs.pctile$tmax90, ">") * ci@namask.mon$tmax) }

## WSDI: Annual.
climdex.wsdi <- function(ci) { return(threshold.exceedance.duration.index(ci@tmax, ci@annual.factor, ci@bs.pctile$tmax90, ">") * ci@namask.mon$tmax) }

## CSDI: Annual.
climdex.csdi <- function(ci) { return(threshold.exceedance.duration.index(ci@tmin, ci@annual.factor, ci@bs.pctile$tmin10, "<") * ci@namask.mon$tmax) }

## DTR: Monthly. Done
climdex.dtr <- function(ci) { return(mean.daily.temp.range(ci@tmax, ci@tmin, ci@monthly.factor) * ci@namask.mon$tavg) }

## Rx1day: Monthly. Should work. Testing?
climdex.rx1day <- function(ci) { return(max.nday.consec.prec(ci@prec, ci@monthly.factor, 1) * ci@namask.mon$prec) }

## Rx5day: Monthly. Should work. Testing?
climdex.rx5day <- function(ci) { return(max.nday.consec.prec(ci@prec, ci@monthly.factor, 5) * ci@namask.mon$prec) }

## SDII: Annual. Should work.
climdex.sdii <- function(ci) { return(simple.precipitation.intensity.index(ci@prec, ci@annual.factor) * ci@namask.ann$prec) }

## R10mm: Annual. Should work.
climdex.r10mm <- function(ci) { return(number.days.op.threshold(ci@prec, ci@annual.factor, 10, ">=") * ci@namask.ann$prec) }

## R20mm: Annual. Should work.
climdex.r20mm <- function(ci) { return(number.days.op.threshold(ci@prec, ci@annual.factor, 20, ">=") * ci@namask.ann$prec) }

## Rnnmm: Annual. Should work.
climdex.rnnmm <- function(ci, threshold) { return(number.days.op.threshold(ci@prec, ci@annual.factor, threshold, ">=") * ci@namask.ann$prec) }

## CDD: Annual. Should work.
climdex.cdd <- function(ci) { return(max.length.spell(ci@prec, ci@annual.factor, 1, "<") * ci@namask.ann$prec) }

## CWD: Annual. Should work.
climdex.cwd <- function(ci) { return(max.length.spell(ci@prec, ci@annual.factor, 1, ">=") * ci@namask.ann$prec) }

## R95pTOT: Annual.
climdex.r95ptot <- function(ci) { return(total.precip.above.threshold(ci@prec, ci@annual.factor, ci@bs.pctile$precwet95) * ci@namask.ann$prec) }

## R99pTOT: Annual.
climdex.r99ptot <- function(ci) { return(total.precip.above.threshold(ci@prec, ci@annual.factor, ci@bs.pctile$precwet99) * ci@namask.ann$prec) }

## PRCPTOT: Annual. Should work.
climdex.prcptot <- function(ci) { return(tapply(ci@prec, ci@annual.factor, sum) * ci@namask.ann$prec) }


##
## HELPERS FINISHED. IMPLEMENTATIION BELOW.
##


## FD, ID, SU, TR, R10mm, R20mm, Rnnmm
number.days.op.threshold <- function(temp, date.factor, threshold, op="<") {
  stopifnot(is.numeric(temp))
  return(tapply(op(temp, threshold), date.factor, sum, na.rm=TRUE))
}

## GSL
## Meaningless if not annual
## Time series must be contiguous
## NOTE: There is a difference of 1 between our output and fclimdex. See line 637; consider case where start and end day are same. Correct answer is 1 day GSL; their answer is 0 day.
growing.season.length <- function(daily.mean.temp, date.factor,
                                  min.length=6, t.thresh=5) {
  return(tapply(daily.mean.temp, date.factor, function(ts) {
    ts.len<- length(ts)
    ts.mid <- floor(ts.len / 2)
    gs.begin <- which(select.blocks.gt.length(ts > t.thresh, min.length - 1))

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
percent.days.op.threshold <- function(temp, date.factor, threshold, op='<') {
  namask <- !is.na(temp)
  return(tapply(op(temp, threshold), date.factor, function(x) { x.nona <- x[!is.na(x)]; if(!length(x.nona)) return(NA); return(sum(x.nona) / length(x.nona) * 100) } ))
}

## WSDI, CSDI
## Thresholds appear to be for each block of 5 days of a year...
threshold.exceedance.duration.index <- function(daily.max.temp, date.factor, warm.thresholds, op=">", min.length=6) {
  periods <- select.blocks.gt.length(op(daily.max.temp, warm.thresholds), min.length - 1)
  return(tapply(periods, date.factor, sum))
}

## DTR
## Max and min temps are assumed to be same length
mean.daily.temp.range <- function(daily.max.temp, daily.min.temp, date.factor) {
  return(tapply(daily.max.temp - daily.min.temp, date.factor, mean))
}

## Rx1day, Rx5day
max.nday.consec.prec <- function(daily.prec, date.factor, ndays) {
  if(ndays == 1) {
    return(tapply(daily.prec, date.factor, max))
  } else {
    ## Ends of the data will be de-emphasized (padded with zero precip data)
    prec.runsum <- runmean(c(rep(0, floor(ndays / 2)), daily.prec, rep(0, floor(ndays / 2))), k=ndays, endrule="trim") * ndays
    return(tapply(prec.runsum, date.factor, max))
  }
}

## SDII
## Period for computation of number of wet days shall be the entire range of the data supplied.
simple.precipitation.intensity.index <- function(daily.prec, date.factor) {
  return(tapply(daily.prec, date.factor, function(prec) { idx <- prec >= 1; return(sum(prec[idx]) / sum(idx)) } ))
}

## CDD, CWD
max.length.spell <- function(daily.prec, date.factor, threshold, op) {
  return(tapply(op(daily.prec, threshold), date.factor, function(x) { return(max(sequential(x))) } ))
}

## R95pTOT, R99pTOT
total.precip.above.threshold <- function(daily.prec, date.factor, threshold) {
  return(tapply(daily.prec[daily.prec > threshold], date.factor, function(x) { return(sum(x, na.rm=TRUE)) } ))
}

## Gotta test this
running.quantile <- function(data, f, n, q, include.mask=NULL) {
  ## Create n lists of indices
  indices.list <- lapply((1:n) - ceiling(n / 2), function(x, indices) { return(indices[max(1, x + 1):min(length(indices), length(indices) + x)]) }, 1:length(data))
  repeated.data <- data[unlist(indices.list)]

  ## Create mask
  bad.mask <- !is.na(repeated.data)
  if(!is.null(include.mask))
    bad.mask <- bad.mask & include.mask[unlist(indices.list)]

  ## Reversing the indices creates the shifted window.
  repeated.f <- f[unlist(rev(indices.list))]
  
  return(t(do.call(data.frame, tapply(repeated.data[bad.mask], repeated.f[bad.mask], quantile, q, type=8))))
}

## Takes a list of booleans; returns a list of booleans where only blocks of TRUE longer than n are still TRUE
select.blocks.gt.length <- function(d, n) {
  if(n == 0)
    return(d)

  if(n >= length(d))
    return(rep(FALSE, length(d)))

  d2 <- Reduce(function(x, y) { return(c(rep(0, y), d[1:(length(d) - y)]) & x) }, 1:n, d)
  return(Reduce(function(x, y) { return(c(d2[(y + 1):length(d2)], rep(0, y)) | x) }, 1:n, d2))
}

## Input vector of booleans
## Returns a vector of integers representing the _length_ of each consecutive sequence of True values
sequential <- function(v) {
  if (! any(v, na.rm=T)) return(0)
  vect <- which(v)
  diff(which(c(T, diff(vect) != 1, T)))
}
