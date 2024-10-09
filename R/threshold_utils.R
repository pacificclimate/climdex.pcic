## Get number of days within range
get.num.days.in.range <- function(x, date.range) {
  return(sum(x >= date.range[1] & x <= date.range[2]))  
}

# Helper function to extract parameters for rxnday indices.
get.rxnday.params <- function(ci, freq= c("monthly", "annual", "seasonal")) {
  stopifnot(!is.null(ci@data$prec))
  data <- ci@data$prec
  date.factors <- ci@date.factors[[match.arg(freq)]]
  mask <- ci@namasks[[match.arg(freq)]]$prec
  cal <- attr(ci@dates, "cal")
  return(list(data = data, date.factor = date.factors, mask = mask, cal = cal))
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
#' ## Calculate frost days.
#' fd <- number.days.op.threshold(ci@@data$tmin,
#'                                ci@@date.factors$annual, 0, "<")
#' 
#' @export
number.days.op.threshold <- function(temp, date.factor, threshold, op="<") {
  stopifnot(is.numeric(temp) && is.numeric(threshold) && is.factor(date.factor))
  return(tapply.fast(match.fun(op)(temp, threshold), date.factor, sum, na.rm=TRUE))
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
#' @param max.missing.days Maximum number of NA values per time period.
#' @return A vector consisting of the mean fraction of days above or below the
#' supplied set of thresholds.
#' @note If date.factor is omitted, daily series will be returned.
#' @seealso \link{climdexInput-class}.
#' @keywords ts climate
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
#' ## Compute monthly tx90p.
#' tx90p <- percent.days.op.threshold(ci@@data$tmax, ci@@dates, ci@@jdays,
#'                                    ci@@date.factors$monthly,
#'                                    ci@@quantiles$tmax$outbase$q90,
#'                                    ci@@quantiles$tmax$inbase$q90,
#'                                    ci@@base.range, ">",
#'                                    ci@@max.missing.days['monthly']) *
#'          ci@@namasks$monthly$tmax
#' 
#' @export
percent.days.op.threshold <- function(temp, dates, jdays, date.factor, threshold.outside.base, base.thresholds, base.range, op='<', max.missing.days) {
  f <- match.fun(op)
  dat <- f(temp, threshold.outside.base[jdays])
  
  inset <- dates >= base.range[1] & dates <= base.range[2]
  ## Don't use in-base thresholds with data shorter than two years; no years to replace with.
  if(sum(inset) > 0 && length(dates) >= 360 * 2) {
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
    dat[inset] <- rowSums(f.result, na.rm=TRUE) / (byrs - 1)
  }
  dat[is.nan(dat)] <- NA
  if(missing(date.factor))
    return(dat)
  na.mask <- get.na.mask(dat, date.factor, max.missing.days)
  ## FIXME: Need to monthly-ize the NA mask calculation, which will be ugly.
  ret <- tapply.fast(dat, date.factor, mean, na.rm=TRUE) * 100 * na.mask
  ret[is.nan(ret)] <- NA
  return(ret)
}

#' @title Sum of spell lengths exceeding daily threshold
#' 
#' @description
#' This function returns the number of spells of more than \code{min.length}
#' days which exceed or are below the given threshold.
#' 
#' @details
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
#' @param max.missing.days Maximum number of NA values per time period.
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
#' phony.date.factor, rep(1:5, 2), rep(1, 5), ">=", 2, TRUE, 1)
#' 
#' ## Without spells spanning years...
#' tedi <- threshold.exceedance.duration.index(prec.dat, phony.date.factor,
#' rep(1:5, 2), rep(1, 5), ">=", 2, FALSE, 1)
#' 
#' @export
threshold.exceedance.duration.index <- function(daily.temp, date.factor, jdays, thresholds, op=">", min.length=6, spells.can.span.years=TRUE, max.missing.days) {
  stopifnot(is.numeric(c(daily.temp, thresholds, min.length)), is.factor(date.factor),
            is.function(match.fun(op)),
            min.length > 0)
  f <- match.fun(op)
  na.mask <- get.na.mask(is.na(daily.temp + thresholds[jdays]), date.factor, max.missing.days)
  
  if(spells.can.span.years) {
    periods <- select.blocks.gt.length(f(daily.temp, thresholds[jdays]), min.length - 1)
    return(tapply.fast(periods, date.factor, sum) * na.mask)
  } else {
    ## fclimdex behaviour...
    return(tapply.fast(1:length(daily.temp), date.factor, function(idx) { sum(select.blocks.gt.length(f(daily.temp[idx], thresholds[jdays[idx]]), min.length - 1)) } ) * na.mask)
  }
}

#' @title Number of days (less than, greater than, etc) a threshold
#' 
#' @description
#' Produces sums of values that exceed (or are below) the specified threshold.
#' 
#' @details
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
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @param mask The mask to be applied to the result.
#' @param freq The frequency of the data (monthly, annual, seasonal).
#' @param cal Calendar object specifying the calendar system.
#' @return A vector containing the the maximum n-day sum of precipitation per
#' time interval. If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.

#' @keywords ts climate
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
#' ## Compute rx5day on a monthly basis.
#' rx5day <- nday.consec.prec.max(ci@@data$prec, ci@@date.factors$monthly, 5)
#' 
#' @export
nday.consec.prec.max <- function(daily.prec, date.factor, ndays, center.mean.on.last.day = FALSE, include.exact.dates = FALSE, mask = 1, freq, cal) {
  stat <- "max"
  if (ndays == 1) {
    if (include.exact.dates) {
      
      df <- exact.date(stat, daily.prec, date.factor, freq, cal, mask)
      return(df)
    }
    return(suppressWarnings(tapply.fast(daily.prec, date.factor, stat, na.rm = TRUE)) * mask)
  }
  ## Ends of the data will be de-emphasized (padded with zero precip data); NAs replaced with 0
  daily.prec[is.na(daily.prec)] <- 0
  prec.runsum <- running.mean(daily.prec, ndays)
  prec.runsum[is.na(prec.runsum)] <- 0
  
  if (center.mean.on.last.day) {
    k2 <- ndays %/% 2
    prec.runsum <- c(rep(0, k2), prec.runsum[1:(length(prec.runsum) - k2)])
  }
  if (include.exact.dates) {
    df <- exact.date(stat, prec.runsum, date.factor, freq, cal, mask)
    df$val <- df$val * ndays
    return(df)
  }
  return((tapply.fast(prec.runsum, date.factor, stat) * ndays) * mask)
}

#' @title Maximum spell length
#' 
#' @description
#' This function returns the longest string of days which exceed or are below
#' the given threshold.
#' 
#' @details
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
#' @param include.exact.dates Logical, if TRUE, return a data frame with spell durations and the start and end dates of each spell; if FALSE, return only the spell durations.
#' @param mask Binary mask indicating valid date.factors.
#' @param cal Calendar object specifying the calendar system.
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

spell.length.max <- function(daily.prec, date.factor, threshold, op, spells.can.span.years, include.exact.dates = FALSE, mask = 1, cal= "365") {
  bools <- match.fun(op)(daily.prec, threshold)
  spells <- get.series.lengths.at.ends(bools)
  if (spells.can.span.years) {
    all.true <- tapply.fast(bools, date.factor, all)
    max.spell <- tapply.fast(spells, date.factor, max)
    
    ## Mask out values which are in the middle of a spell with NA
    na.mask <- c(1, NA)[as.integer((max.spell == 0) & all.true) + 1]
    max.spell <- max.spell * na.mask
  } else {
    max.spell <- tapply.fast(bools, date.factor, function(x) {
      max(get.series.lengths.at.ends(x))
    })
  }
  if (include.exact.dates) {
    end <- tapply.fast(bools, date.factor, function(x) {
      which.max(get.series.lengths.at.ends(x))
    })
    
    start <- end - max.spell
    origin <- paste(names(max.spell), "01", "01", sep = "-")
    origin.pcict <- as.PCICt(origin, cal)
    seconds.per.day <- 86400
    start.pcict <- origin.pcict + start * seconds.per.day
    end.pcict <- origin.pcict + (end - 1) * seconds.per.day
    
    df <- data.frame(
      start = format(start.pcict, "%Y-%m-%d"),
      duration = max.spell * mask,
      end = format(end.pcict, "%Y-%m-%d")
    )
    
    df$start[which(df$duration == 0)]<- NA
    df$end[which(df$duration == 0)] <- NA
    
    df$start[is.na(df$duration)] <- NA
    df$end[is.na(df$duration)] <- NA
    return(df)
  } else {
    return(max.spell * mask)
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