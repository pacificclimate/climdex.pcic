#' Any of tmin (daily minimum temperature), tmax (daily maximum temperature),
#' tavg (daily mean temperature), and prec (daily precipitation) can be passed
#' in. tavg will be derived from the mean of tmax and tmin if it is not
#' supplied. If any of tmin, tmax, and prec are not supplied, the set of
#' indices which can be calculated will be limited to indices which do not
#' involve the missing variables.
#'
#' For all data supplied, the associated dates must also be supplied.
#' 
