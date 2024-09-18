#' @importFrom stats quantile
#' 

## Check that climdexInput data structure is valid.
valid.climdexInput <- function(object) {
  temp.quantiles <- c(10, 90)
  prec.quantiles <- c(95, 99)
  errors <- c()
  
  separate.base <- c(tmax=T, tmin=T, tavg=T, prec=F)
  present.data.vars <- names(object@data)
  length.check.slots <- c("dates", "jdays")
  length.check.members <- c("date.factors", "data")
  data.lengths <- c(sapply(object@data, length), sapply(length.check.slots, function(y) length(slot(object, y))), unlist(sapply(length.check.members, function(y) { sapply(slot(object, y), length) })))
  quantiles <- list(tmax=temp.quantiles, tmin=temp.quantiles, prec=prec.quantiles)
  
  if(!all(data.lengths == max(data.lengths)))
    errors <- c(errors, "Data fields, dates, and date factors must all be of the same length")
  
  ## Check that namasks have columns for each of the variables
  if(!all(c("annual", "monthly", "seasonal") %in% names(object@namasks)) || !all(present.data.vars %in% names(object@namasks$annual) & present.data.vars %in% names(object@namasks$monthly)& present.data.vars %in% names(object@namasks$seasonal)))
    errors <- c(errors, "NA mask for monthly, seasonal and annual must contain data for all variables supplied.")
  
  ## Check that appropriate thresholds are present.
  need.base.data <- get.num.days.in.range(object@dates, object@base.range) > 0
  
  if(length(object@northern.hemisphere) != 1)
    errors <- c(errors, "northern.hemisphere must be of length 1.")
  
  if(length(errors) == 0)
    return(TRUE)
  else
    return(errors)
}

## Returns PCICt field or dies
get.date.field <- function(input.data, cal, date.types) {
  valid.date.types <- sapply(date.types, function(x) { return(!inherits(try(input.data[,x$fields], silent=TRUE), "try-error")) })
  
  if(sum(valid.date.types) == 0) {
    stop("Could not find a workable set of date fields")
  }
  
  date.type <- date.types[[which(valid.date.types)[1]]]
  date.strings <- do.call(paste, input.data[,date.type$fields,drop=FALSE])
  return(as.PCICt(date.strings, format=date.type$format, cal=cal))
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
  
  if(all(c(is.null(tmax), is.null(tmin), is.null(prec), is.null(tavg))))
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
  
  if(!inherits(quantiles, "list"))
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


#' Get available indices by name
#'
#' This function returns a vector of (function) names of available indices.
#'
#' This function takes a climdexInput object as input and returns the names of
#' all the indices which may be computed or, if \code{get.function.names} is
#' TRUE (the default), the names of the functions corresponding to the indices.
#' 
#' @param ci Object of type climdexInput.
#' @param function.names Whether to return function names.
#' @return A vector containing an annual timeseries of precipitation in wet days.
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
#' ## Get list of functions which might be run.
#' func.names <- climdex.get.available.indices(ci)
#' @export
climdex.get.available.indices <- function(ci, function.names=TRUE) {
  available.indices <- list(tmax=c('su', 'id', 'txx', 'txn', 'tx10p', 'tx90p', 'wsdi'),
                            tmin=c('fd', 'tr', 'tnx', 'tnn', 'tn10p', 'tn90p', 'csdi'),
                            tavg=c('gsl', 'dtr'),
                            prec=c('rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot'))
  if(function.names) {
    return(paste("climdex", unlist(available.indices[names(ci@data)]), sep="."))
  } else {
    return(unlist(available.indices[names(ci@data)], use.names=FALSE))
  }
}
