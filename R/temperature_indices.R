#' Monthly Maximum of Daily Maximum Temperature
#'
#' This function computes the climdex index TXx.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly, seasonal or annual maximum of daily maximum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @return A vector containing the value of the index per time interval. 
#' If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.
#' @template generic_seealso_references
#' @templateVar cdxvar txx
#' @templateVar cdxdescription a monthly timeseries of maximum daily maximum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.txx <- function(ci, freq = c("monthly", "annual", "seasonal"), include.exact.dates = FALSE) {
  return(compute.stat(ci, "max", "tmax", freq, include.exact.dates))
}

#' Monthly Maximum of Daily Minimum Temperature
#'
#' This function computes the climdex index TNx.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly, seasonal or annual maximum of daily minimum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @return A vector containing the value of the index per time interval. 
#' If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.
#' @template generic_seealso_references
#' @templateVar cdxvar tnx
#' @templateVar cdxdescription a monthly timeseries of maximum daily minimum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.tnx <- function(ci, freq = c("monthly", "annual", "seasonal"), include.exact.dates = FALSE) {
  return(compute.stat(ci, "max", "tmin", freq, include.exact.dates))
}

#' Monthly Minimum of Daily Maximum Temperature
#'
#' This function computes the climdex index TXn.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly, seasonal or annual minimum of daily maximum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @return A vector containing the value of the index per time interval. 
#' If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.
#' @template generic_seealso_references
#' @templateVar cdxvar txn
#' @templateVar cdxdescription a monthly timeseries of minimum daily maximum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.txn <- function(ci, freq = c("monthly", "annual", "seasonal"), include.exact.dates = FALSE) {
  return(compute.stat(ci, "min", "tmax", freq, include.exact.dates))
}

#' Monthly Minimum of Daily Minimum Temperature
#'
#' This function computes the climdex index TNn.
#' 
#' This function takes a climdexInput object as input and computes
#' the monthly, seasonal or annual minimum of daily minimum temperature.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the index values and exact dates for each month, season, or year; if FALSE, return only the index values.
#' @return A vector containing the value of the index per time interval. 
#' If `include.exact.dates` is TRUE, a data frame with additional exact dates is returned.
#' @template generic_seealso_references
#' @templateVar cdxvar tnn
#' @templateVar cdxdescription a monthly timeseries of minimum daily minimum temperature.
#' @template get_generic_example
#' 
#' @export
climdex.tnn <- function(ci, freq = c("monthly", "annual", "seasonal"), include.exact.dates = FALSE) {
  return(compute.stat(ci, "min", "tmin", freq, include.exact.dates))
}

## Our implementation currently follows the example set by fclimdex for dealing with missing values, which is wrong; it biases results upwards when missing values are present.

#' Percent of Values Below 10th Percentile Daily Minimum Temperature
#' 
#' This function computes the climdex index TN10p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly, seasonal or annual percent of values below the 10th percentile of baseline
#' daily minimum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tn10p
#' @templateVar cdxdescription a monthly timeseries of the TN10p index.
#' @template get_generic_example
#' 
#' @export
climdex.tn10p <- function(ci, freq=c("monthly", "annual", "seasonal")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(percent.days.op.threshold(ci@data$tmin, ci@dates, ci@jdays, ci@date.factors[[match.arg(freq)]], ci@quantiles$tmin$outbase$q10, ci@quantiles$tmin$inbase$q10, ci@base.range, "<", ci@max.missing.days[match.arg(freq)]) * ci@namasks[[match.arg(freq)]]$tmin) }

#' Percent of Values Below 10th Percentile Daily Maximum Temperature
#' 
#' This function computes the climdex index TX10p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly, seasonal or annual percent of values below the 10th percentile of baseline
#' daily maximum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tx10p
#' @templateVar cdxdescription a monthly timeseries of the TX10p index.
#' @template get_generic_example
#' 
#' @export
climdex.tx10p <- function(ci, freq=c("monthly", "annual", "seasonal")) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, ci@date.factors[[match.arg(freq)]], ci@quantiles$tmax$outbase$q10, ci@quantiles$tmax$inbase$q10, ci@base.range, "<", ci@max.missing.days[match.arg(freq)]) * ci@namasks[[match.arg(freq)]]$tmax) }

#' Percent of Values Above 90th Percentile Daily Minimum Temperature
#' 
#' This function computes the climdex index TN90p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly, seasonal or annual percent of values above the 90th percentile of baseline
#' daily minimum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tn90p
#' @templateVar cdxdescription a monthly timeseries of the TN90p index.
#' @template get_generic_example
#' 
#' @export
climdex.tn90p <- function(ci, freq=c("monthly", "annual", "seasonal")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(percent.days.op.threshold(ci@data$tmin, ci@dates, ci@jdays, ci@date.factors[[match.arg(freq)]], ci@quantiles$tmin$outbase$q90, ci@quantiles$tmin$inbase$q90, ci@base.range, ">", ci@max.missing.days[match.arg(freq)]) * ci@namasks[[match.arg(freq)]]$tmin) }

#' Percent of Values Above 90th Percentile Daily Maximum Temperature
#' 
#' This function computes the climdex index TX90p.
#' 
#' This function takes a climdexInput object as input and computes the
#' monthly, seasonal or annual percent of values above the 90th percentile of baseline
#' daily maximum temperature.
#' 
#' @template threshold_indices_block
#' @template threshold_indices_args
#' @template missing_values_caveat
#'
#' @template generic_seealso_references
#' @templateVar cdxvar tx90p
#' @templateVar cdxdescription a monthly timeseries of the TX90p index.
#' @template get_generic_example
#' 
#' @export
climdex.tx90p <- function(ci, freq=c("monthly", "annual", "seasonal")) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(percent.days.op.threshold(ci@data$tmax, ci@dates, ci@jdays, ci@date.factors[[match.arg(freq)]], ci@quantiles$tmax$outbase$q90, ci@quantiles$tmax$inbase$q90, ci@base.range, ">", ci@max.missing.days[match.arg(freq)]) * ci@namasks[[match.arg(freq)]]$tmax) }

#' Mean Diurnal Temperature Range
#' 
#' This function computes the diurnal temperature range on a monthly basis.
#' 
#' \code{climdex.dtr} computes the mean daily diurnal temperature range. The
#' frequency of observation can be either monthly or annual.
#' 
#' @param ci Object of type climdexInput.
#' @param freq Time frequency to aggregate to.
#' @return A vector containing the mean monthly, mean seasonal or mean annual diurnal
#' temperature range.
#' @note This function creates results which may differ in the 3rd decimal
#' place from the results from fclimdex.
#' @template generic_seealso_references
#' @templateVar cdxvar dtr
#' @templateVar cdxdescription a monthly timeseries of mean diurnal temperature range.
#' @template get_generic_example
#' 
#' @export
climdex.dtr <- function(ci, freq=c("monthly", "annual", "seasonal")) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@data$tmax) && !is.null(ci@data$tavg)); return(compute.mean.daily.temp.range(ci@data$tmax, ci@data$tmin, ci@date.factors[[match.arg(freq)]]) * ci@namasks[[match.arg(freq)]]$tavg) }

## DTR
## Computes mean diurnal temperature range in each period (as specified by date.factor).
## Max and min temps are assumed to be same length
compute.mean.daily.temp.range <- function(daily.max.temp, daily.min.temp, date.factor) {
  dat <- tapply.fast(daily.max.temp - daily.min.temp, date.factor, mean, na.rm=TRUE)
  dat[is.nan(dat)] <- NA
  dat
}