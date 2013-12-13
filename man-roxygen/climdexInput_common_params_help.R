#' The \code{base.range} argument is a pair of 4 digit years which bound the
#' data on which the base percentiles are calculated.
#' 
#' The \code{n} argument specifies the size of the window used when computing
#' the percentiles used in \code{\link{climdex.tx10p}},
#' \code{\link{climdex.tn10p}}, \code{\link{climdex.tx90p}}, and
#' \code{\link{climdex.tn90p}}.
#' 
#' The \code{northern.hemisphere} argument specifies whether the data came from
#' the northern hemisphere. If FALSE, data is assumed to have come from the
#' southern hemisphere. This is used when computing growing season length; if
#' the data is from the southern hemisphere, growing season length is the
#' growing season starting in the beginning of July of the year indicated,
#' running to the end of June of the following year.
#' 
#' The \code{pad.data.with.first.last.values} argument specifies whether to pad
#' the data passed into the baseline quantile routine with the first and last
#' values. If TRUE, the first (at the beginning of the series) and last (at the
#' end of the series) values will be used to pad the beginning and ends of this
#' series. If FALSE, either NA or the values for the previous two (at the
#' beginning) and last two (at the end) days of data will be used.
#' 
#' The \code{quantiles} argument allows the user to supply pre-computed quantiles.
#' This is a list consisting of quantiles for each variable.
#' 
#' For each temperature variable, there are separate lists of quantiles for 
#' inbase and outbase, with these names. In both cases, quantiles within these
#' lists are named q10 for the 10th percentile and q90 for the 90th percentile.
#' Other percentiles would be named qnn for the nnth percentile. For the
#' outbase quantiles, each element in the list is a vector of length 365 (or 360
#' in the case of 360-day calendars), corresponding to one value for each day of
#' the year. For the inbase quantiles, each element in the list is an array of
#' dimensions [365 or 360, nyr, nyr - 1], where nyr is the number of years in
#' the base period. Each value corresponds to a quantile for each day, for each
#' year, with a particular year replaced.
#'
#' For precipitation variables, there is a named vector of quantiles, consisting
#' of at least q95 and q99. 
#'
#' The \code{temp.qtiles} and \code{prec.qtiles} arguments allow the user to
#' modify the quantiles calculated. For example, specifying
#' temp.qtiles=c(0.10, 0.50, 0.90) would calculate the 10th, 50th, and 90th
#' percentiles for temperature.
#' 
#' @seealso \code{\link{climdex.pcic-package}}, \code{\link{strptime}}.
#' @references \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
#' @keywords ts climate
