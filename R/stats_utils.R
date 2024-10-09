## Lower overhead version of tapply
tapply.fast <- function(X, INDEX, FUN = NULL, ..., simplify = TRUE) {
  FUN <- if (!is.null(FUN)) {
    match.fun(FUN)
  }
  
  if (!is.factor(INDEX)) {
    stop("INDEX must be a factor.")
  }
  
  if (length(INDEX) != length(X)) {
    stop("arguments must have the same length")
  }
  
  if (is.null(FUN)) {
    return(INDEX)
  }
  
  
  ans <- lapply(split(X, INDEX), FUN, ...)
  if (is.function(FUN) && (identical(FUN, which.min) || identical(FUN, which.max))) {
    # Handle which.min & which.max separately
    ans <- lapply(ans, function(x) if (length(x) == 0) NA else x)
    ans <- unlist(ans)
  } else {
    ans <- unlist(ans, recursive = FALSE)
  }
  names(ans) <- levels(INDEX)
  
  return(ans)
}

## Calculate a running quantile on the data set over the bootstrap range.
## If get.bootstrap.data is TRUE, use the Zhang boostrapping method described in Xuebin Zhang et al's 2005 paper, "Avoiding Inhomogeneity in Percentile-Based Indices of Temperature Extremes" J.Clim vol 18 pp.1647-1648, "Removing the 'jump'".
## Expects PCICt for all dates
zhang.running.qtile <- function(x, dates.base, qtiles, bootstrap.range, include.mask=NULL, n=5, get.bootstrap.data=FALSE, min.fraction.present=0.1) {
  inset <- get.bootstrap.set(dates.base, bootstrap.range, n)
  dpy <- ifelse(is.null(attr(dates.base, "dpy")), 365, attr(dates.base, "dpy"))
  nyears <- floor(sum(inset) / dpy)
  
  if(!is.null(include.mask))
    x[include.mask] <- NA
  
  bs.data <- x[inset]
  
  qdat <- NULL
  if(get.bootstrap.data) {
    d <- .Call("running_quantile_windowed_bootstrap", bs.data, n, qtiles, dpy, min.fraction.present, PACKAGE='climdex.pcic')
    dim(d) <- c(dpy, nyears, nyears - 1, length(qtiles))
    qdat <- lapply(1:length(qtiles), function(x) { r <- d[,,,x, drop=FALSE]; dim(r) <- dim(r)[1:3]; r })
  } else {
    res <- running.quantile(bs.data, n, qtiles, dpy, min.fraction.present)
    qdat <- lapply(1:length(qtiles), function(x) { res[,x] })
  }
  names(qdat) <- paste("q", qtiles * 100, sep="")
  return(qdat)
}

get.temp.var.quantiles <- function(filled.data, date.series, bs.date.series, qtiles, bs.date.range, n, in.base=FALSE, min.base.data.fraction.present=0.1) {
  base.data <- create.filled.series(filled.data, date.series, bs.date.series)
  if(in.base)
    return(list(outbase=zhang.running.qtile(base.data, dates.base=bs.date.series, qtiles=qtiles, bootstrap.range=bs.date.range, n=n, min.fraction.present=min.base.data.fraction.present),
                inbase=zhang.running.qtile(base.data, dates.base=bs.date.series, qtiles=qtiles, bootstrap.range=bs.date.range, n=n, get.bootstrap.data=TRUE, min.fraction.present=min.base.data.fraction.present)))
  else
    return(list(outbase=zhang.running.qtile(base.data, dates.base=bs.date.series, qtiles=qtiles, bootstrap.range=bs.date.range, n=n, min.fraction.present=min.base.data.fraction.present)))
}

get.prec.var.quantiles <- function(filled.prec, date.series, bs.date.range, qtiles=c(0.95, 0.99)) {
  wet.days <- !(is.na(filled.prec) | filled.prec < 1)
  inset <- date.series >= bs.date.range[1] & date.series <= bs.date.range[2] & !is.na(filled.prec) & wet.days
  pq <- quantile(filled.prec[inset], qtiles, type=8)
  names(pq) <- paste("q", qtiles * 100, sep="")
  return(pq)
}

#' @title Method for getting threshold quantiles for use in computing indices
#' 
#' @description
#' This function creates threshold quantiles for use with climdexInput.raw
#' or climdexInput.csv.
#' 
#' @details
#' This function takes input climate data at daily resolution, and produces as
#' output a set of threshold quantiles. This data structure can then be passed
#' to climdexInput.raw or climdexInput.csv.
#'
#' For more details on arguments, see \code{\link{climdexInput.raw}}.
#'
#' @seealso \code{\link{climdex.pcic-package}}, \code{\link{climdexInput.raw}}.
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @keywords ts climate
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
#' quantiles <- get.outofbase.quantiles(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#'
#' @export
get.outofbase.quantiles <- function(tmax=NULL, tmin=NULL, prec=NULL, tmax.dates=NULL, tmin.dates=NULL, prec.dates=NULL, base.range=c(1961, 1990), n=5, temp.qtiles=c(0.10, 0.90), prec.qtiles=c(0.95, 0.99), min.base.data.fraction.present=0.1) {
  days.threshold <- 359
  check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n)
  
  d.list <- list(tmin.dates, tmax.dates, prec.dates)
  all.dates <- do.call(c, d.list[!sapply(d.list, is.null)])
  last.day.of.year <- get.last.monthday.of.year(all.dates)
  cal <- attr(all.dates, "cal")
  
  bs.date.range <- as.PCICt(paste(base.range, c("01-01", last.day.of.year), sep="-"), cal=cal)
  new.date.range <- as.PCICt(paste(as.numeric(format(range(all.dates), "%Y", tz="GMT")), c("01-01", last.day.of.year), sep="-"), cal=cal)
  date.series <- seq(new.date.range[1], new.date.range[2], by="day")
  bs.date.series <- seq(bs.date.range[1], bs.date.range[2], by="day")
  
  quantiles <- list()
  
  if(!is.null(tmax)) {
    if(get.num.days.in.range(tmax.dates, bs.date.range) <= days.threshold)
      stop("There is less than a year of tmax data within the base period. Consider revising your base range and/or check your input data.")
    filled.tmax <- create.filled.series(tmax, trunc(tmax.dates, "days"), date.series)
    quantiles$tmax <- get.temp.var.quantiles(filled.tmax, date.series, bs.date.series, temp.qtiles, bs.date.range, n)
  } 
  
  if(!is.null(tmin)) {
    if(get.num.days.in.range(tmin.dates, bs.date.range) <= days.threshold)
      stop("There is less than a year of tmin data within the base period. Consider revising your base range and/or check your input data.")
    filled.tmin <- create.filled.series(tmin, trunc(tmin.dates, "days"), date.series)
    quantiles$tmin <- get.temp.var.quantiles(filled.tmin, date.series, bs.date.series, temp.qtiles, bs.date.range, n)
  }
  
  if(!is.null(prec)) {
    if(get.num.days.in.range(prec.dates, bs.date.range) <= days.threshold)
      stop("There is less than a year of prec data within the base period. Consider revising your base range and/or check your input data.")
    filled.prec <- create.filled.series(prec, trunc(prec.dates, "days"), date.series)
    quantiles$prec <- get.prec.var.quantiles(filled.prec, date.series, bs.date.range, prec.qtiles)
  }
  return(quantiles)
}

# Computes a specified statistic (min, max) for a given climate index and frequency.
compute.stat <- function(ci, stat, data.key, freq = c("monthly", "annual", "seasonal"), include.exact.dates) {
  stopifnot(!is.null(ci@data[[data.key]]))
  data <- ci@data[[data.key]]
  date.factors <- ci@date.factors[[match.arg(freq)]]
  mask <- ci@namasks[[match.arg(freq)]][[data.key]]
  cal <- attr(ci@dates, "cal")
  
  if (include.exact.dates) {
    return(exact.date(stat, data, date.factors, freq, cal, mask))
  }
  
  return(suppressWarnings(tapply.fast(data, date.factors, stat, na.rm = TRUE)) * mask)
}

## Returns an n-day running quantile for each day of data (dimensions c(dpy, q))
running.quantile <- function(data, n, q, dpy, min.fraction) {
  ret <- .Call("running_quantile_windowed", data, n, q, dpy, min.fraction, PACKAGE='climdex.pcic')
  dim(ret) <- c(length(q), dpy)
  return(t(ret))
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

#' @title Running Mean of a Vector
#'
#' @description Calculates the running means of a vector with a shifting window
#'
#' @details Returns a new vector the same length as vec, where the ith
#' element is the mean of the bin of elements centered at the ith element
#' of the original vector. Means cannot be calculated for elements less
#' than half the width of the bin from the beginning or end of the vector;
#' the result vector has NA in those positions.
#'
#' @param vec A vector
#' @param bin The number of entries to average over for each mean
#'
#' @return a vector containing the running mean of bin elements of vec
#'
#' @examples
#' 
#' \dontrun{
#' running.mean(c(1, 2, 3, 4, 5, 6), 2)
#' }
#' \dontrun{
#' running.mean(c(5, 5, 5, 5, 5), 4)
#' }
running.mean <- function(vec, bin){
  vec = as.vector(vec)
  len = length(vec)
  if (bin<=1) {
    return (vec)
  }
  if (bin > len) {
    bin = len
  }
  left.bin = bin%/%2
  
  means = double(len)
  
  right.bin = bin - left.bin - 1
  means = c( sum(vec[1:bin]), diff(vec,bin) ) # find the first sum and the differences from it
  means = cumsum(means)/bin                  # apply precomputed differences
  means = c(rep(NA,left.bin), means, rep(NA,right.bin))   # extend to original vector length
  return(means)
}
