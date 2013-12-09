library(climdex.pcic)
library(RUnit)

climdex.pcic.test.all.indices <- function() {
  ## Parse the dates into PCICt.
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  
  ## Load the data in.
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION, tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))

  all.indices <- c('fd', 'su', 'id', 'tr', 'gsl', 'txx', 'tnx', 'txn', 'tnn', 'tn10p', 'tx10p', 'tn90p', 'tx90p', 'wsdi', 'csdi',
                   'dtr', 'rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')

  for(index in all.indices) {
    fun <- match.fun(paste('climdex', index, sep="."))
    valid.result <- get(paste('ec.1018935', index, sep="."))
    checkEquals(valid.result, fun(ci), paste(index, "didn't match expected result"))
  }
}
