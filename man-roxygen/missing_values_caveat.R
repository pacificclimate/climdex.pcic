#' @note These functions may calculate slightly different results than
#' fclimdex.
#' 
#' The bootstrapping method is not well defined for cases where the base data
#' contains numerous missing values.  Because of that, this code (and fclimdex)
#' are not very robust against missing values with respect to these indicies.
#' When computing percentiles inside the base period, both this implementation
#' and fclimdex do not divide through by the number of non-missing values when
#' aggregating the values inside the base period. Instead, they divide through
#' by the number of base years minus one. This will result in a negative bias
#' when missing values are present.
