#' Monthly Maximum 1-day Precipitation
#' 
#' This function computes the climdex index Rx1day.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index Rx1day: monthly, seasonal or annual maximum 1-day precipitation.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @return A vector containing the value of the index per time interval. 
#' If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.
#' @template rx5day_common
#' @template generic_seealso_references
#' @templateVar cdxvar rx1day
#' @templateVar cdxdescription a timeseries of monthly maximum 1-day precipitation.
#' @template get_generic_example
#' 
#' @export
climdex.rx1day <- function(ci, freq = c("monthly", "annual", "seasonal"), include.exact.dates = FALSE) {
  rxnday.params <- get.rxnday.params(ci, freq)
  ndays <- 1
  return(nday.consec.prec.max(daily.prec = rxnday.params$data, date.factor = rxnday.params$date.factor, ndays = ndays, include.exact.dates = include.exact.dates,  mask = rxnday.params$mask, freq = freq, cal = rxnday.params$cal))
}

#' Monthly Maximum 5-day Consecutive Precipitation
#' 
#' This function computes the climdex index Rx5day.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index Rx5day: monthly, seasonal or annual maximum 5-day consecutive precipitation.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param center.mean.on.last.day Whether to center the 5-day running mean on
#' the last day of the window, instead of the center day.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @return A vector containing the value of the index per time interval. 
#' If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.
#' @template rx5day_common
#' @template generic_seealso_references
#' @templateVar cdxvar rx5day
#' @templateVar cdxdescription a timeseries of monthly maximum 5-day consecutive precipitation.
#' @template get_generic_example
#' 
#' @export
climdex.rx5day <- function(ci, freq = c("monthly", "annual", "seasonal"), center.mean.on.last.day = FALSE, include.exact.dates = FALSE) {
  rxnday.params <- get.rxnday.params(ci, freq)
  ndays <- 5
  return(nday.consec.prec.max(daily.prec = rxnday.params$data, date.factor = rxnday.params$date.factor, ndays = ndays, center.mean.on.last.day = center.mean.on.last.day, include.exact.dates = include.exact.dates, mask = rxnday.params$mask, freq = freq, cal = rxnday.params$cal))
}

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
climdex.sdii <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(simple.precipitation.intensity.index(ci@data$prec, ci@date.factors$annual) * ci@namasks$annual$prec) }

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
#' @templateVar cdxdescription an annual timeseries of the R10mm index.
#' @template get_generic_example
#' 
#' @export
climdex.r10mm <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(number.days.op.threshold(ci@data$prec, ci@date.factors$annual, 10, ">=") * ci@namasks$annual$prec) }

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
#' @templateVar cdxdescription an annual timeseries of the R20mm index.
#' @template get_generic_example
#' 
#' @export
climdex.r20mm <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(number.days.op.threshold(ci@data$prec, ci@date.factors$annual, 20, ">=") * ci@namasks$annual$prec) }

#' Precipitation Exceeding A Specified Amount Per Day
#' 
#' This function computes the climdex index Rnnmm.
#'
#' This function takes a climdexInput object as input and computes the climdex
#' index Rnnmm: the annual count of days where daily precipitation is more than \code{nn} mm per day.
#' 
#' @param ci Object of type climdexInput.
#' @param threshold The threshold to be used for Rnnmm.
#' @return A vector containing the value of the index for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar rnnmm
#' @templateVar cdxdescription an annual timeseries of the R1mm index.
#' @template get_generic_example
#' 
#' @export
climdex.rnnmm <- function(ci, threshold=1) {
  stopifnot(!is.null(ci@data$prec));
  if(!is.numeric(threshold) || length(threshold) != 1) stop("Please specify a single numeric threshold value.");
  
  return(number.days.op.threshold(ci@data$prec, ci@date.factors$annual, threshold, ">=") * ci@namasks$annual$prec)
}

#' Maximum Consecutive Dry Days
#' 
#' This function computes the climdex index CDD.
#'
#' This function computes the climdex index CDD: the annual maximum length of dry spells, in days.
#' Dry spells are considered to be sequences of days where daily preciptation
#' is less than 1mm per day.
#' @param include.exact.dates Logical, if TRUE, return a data frame with spell durations and the start and end dates of each spell; if FALSE, return only the spell durations.
#' @template cdd_common
#' @templateVar cdxvar cdd
#' @templateVar cdxdescription an annual timeseries of the CDD index.
#' @template get_generic_example
#' 
#' @export
climdex.cdd <- function(ci, spells.can.span.years = T, include.exact.dates = FALSE) {
  stopifnot(!is.null(ci@data$prec))
  return(spell.length.max(ci@data$prec, ci@date.factors$annual, 1, "<", spells.can.span.years, include.exact.dates, ci@namasks$annual$prec, attr(ci@dates, "cal")))
}

#' Maximum Consecutive Wet Days
#' 
#' This function computes the climdex index CWD.
#' 
#' This function takes a climdexInput object as input and computes the climdex
#' index CWD: the annual maximum length of wet spells, in days.
#' Wet spells are considered to be sequences of days where daily precipitation
#' is at least 1mm per day.
#' @param include.exact.dates Logical, if TRUE, return a data frame with spell durations and the start and end dates of each spell; if FALSE, return only the spell durations.
#' @template cdd_common
#' @templateVar cdxvar cdd
#' @templateVar cdxdescription an annual timeseries of the CWD index.
#' @template get_generic_example
#' 
#' @export
climdex.cwd <- function(ci, spells.can.span.years = T, include.exact.dates = FALSE) {
  stopifnot(!is.null(ci@data$prec))
  return(spell.length.max(ci@data$prec, ci@date.factors$annual, 1, ">=", spells.can.span.years, include.exact.dates, ci@namasks$annual$prec, attr(ci@dates, "cal")))
}

#' Total Daily Precipitation Exceeding 95\\%ile Threshold
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
#' @templateVar cdxdescription an annual timeseries of the R95pTOT index.
#' @template get_generic_example
#' 
#' @export
climdex.r95ptot <- function(ci) { stopifnot(!is.null(ci@data$prec) && !is.null(ci@quantiles$prec)); return(total.precip.op.threshold(ci@data$prec, ci@date.factors$annual, ci@quantiles$prec['q95'], ">") * ci@namasks$annual$prec) }

#' Total Daily Precipitation Exceeding 99\\%ile Threshold
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
#' @templateVar cdxdescription an annual timeseries of the R99pTOT index.
#' @template get_generic_example
#' 
#' @export
climdex.r99ptot <- function(ci) { stopifnot(!is.null(ci@data$prec) && !is.null(ci@quantiles$prec)); return(total.precip.op.threshold(ci@data$prec, ci@date.factors$annual, ci@quantiles$prec['q99'], ">") * ci@namasks$annual$prec) }

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
climdex.prcptot <- function(ci) { stopifnot(!is.null(ci@data$prec)); return(total.precip.op.threshold(ci@data$prec, ci@date.factors$annual, 1, ">=") * ci@namasks$annual$prec) }


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
