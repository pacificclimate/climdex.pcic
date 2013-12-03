#' Computation of these percentiles involves use of a boostrap procedure,
#' described below but described in more depth in [Zhang, 2005].
#' 
#' Computation of these values outside of the base period involves comparing
#' the temperature data for each day with the corresponding percentiles for a 5
#' day running window surrounding that day. The resulting monthly series is
#' then the monthly percentage of values that meet the criteria.
#' 
#' Computation of these values inside the base period is more complicated. It
#' involves comparison of the daily temperature data with the corresponding day
#' of temperature data in each of (n - 1) sets of data. The sets consist of the
#' data for the base period with the current year replaced with each of the
#' other years. The results of these comparisons are then averaged to give a
#' value between 0 and 1. Finally, the resulting daily series is aggregated to
#' a monthly series by averaging these daily values and multiplying by 100 to
#' give a monthly percentile value.
