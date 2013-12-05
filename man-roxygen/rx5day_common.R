#' @return A vector containing the value of the index for each month of each
#' year.
#' @note The default behaviour of climdex.rx5day differs somewhat from
#' fclimdex, as fclimdex and climdex.pcic differ on the definition of Rx5day.
#' The running sum series computed by fclimdex is off by 2 days, and the first
#' day a running sum can be computed for is left out entirely. The behaviour of
#' fclimdex can be replicated by setting center.mean.on.last.day to TRUE.
