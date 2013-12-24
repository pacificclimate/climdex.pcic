#' @param ci Object of type climdexInput.
#' @param spells.can.span.years Whether to allow spells of dry/wet days to span
#' years.
#' @return A vector containing the value of the index for each year.
#' @note These functions may calculate slightly different results than
#' fclimdex.
#' 
#' Behaviour of climdex.wsdi and climdex.csdi differ somewhat from fclimdex.
#' fclimdex considers all days in a warm or cold spell to be part of the year
#' in which the spell ended.  climdex.wsdi and climdex.csdi split the spell
#' such that days in each spell are allocated to the separate years in the days
#' occurred.
#'
#' @seealso \code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}},
#' \code{\link{threshold.exceedance.duration.index}}.
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @keywords ts climate
