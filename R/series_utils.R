## Creates a filled series given the data, dates, and new date sequence to be used.
create.filled.series <- function(data, data.dates, new.date.sequence) {
  new.data <- rep(NA, length(new.date.sequence))
  data.in.new.data <- (data.dates >= new.date.sequence[1]) & (data.dates <= new.date.sequence[length(new.date.sequence)])
  indices <- floor(as.numeric(data.dates[data.in.new.data] - new.date.sequence[1], units="days")) + 1
  new.data[indices] <- data[data.in.new.data]
  return(new.data)
}

## Get set of days for bootstrap use
get.bootstrap.set <- function(dates, bootstrap.range, win.size) {
  dpy <- ifelse(is.null(attr(dates, "dpy")), 365, attr(dates, "dpy"))
  return(dates >= bootstrap.range[1] & dates <= bootstrap.range[2] & (dpy == 360 | format(dates, format="%m-%d", tz="GMT") != "02-29"))
}

## Get NA mask given threshold and split factor
get.na.mask <- function(x, f, threshold) {
  return(c(1, NA)[1 + as.numeric(tapply.fast(is.na(x), f, function(y) { return(sum(y) > threshold) } ))])
}

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
#' series.lengths <- get.series.lengths.at.ends(c(TRUE, TRUE, TRUE, FALSE,
#' TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE))
#' 
#' 
#' @export
get.series.lengths.at.ends <- function(x, na.value=FALSE) {
  stopifnot(is.logical(x) && is.logical(na.value))
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
