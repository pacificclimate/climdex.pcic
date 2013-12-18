library(climdex.pcic)
library(RUnit)

climdex.pcic.test.intake.routines <- function() {
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmax.dat <- ec.1018935.tmax$MAX_TEMP
  tmin.dat <- ec.1018935.tmin$MIN_TEMP
  prec.dat <- ec.1018935.prec$ONE_DAY_PRECIPITATION
  all.indices <- c('fd', 'su', 'id', 'tr', 'gsl', 'txx', 'tnx', 'txn', 'tnn', 'tn10p', 'tx10p', 'tn90p', 'tx90p', 'wsdi', 'csdi',
                   'dtr', 'rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')
  tmax.indices <- c('su', 'id', 'txx', 'txn', 'tx10p', 'tx90p', 'wsdi')
  tmin.indices <- c('fd', 'tr', 'tnx', 'tnn', 'tn10p', 'tn90p', 'csdi')
  tmax.tmin.indices <- c('gsl', 'dtr')
  prec.indices <- c('rx1day', 'rx5day', 'sdii', 'r10mm', 'r20mm', 'rnnmm', 'cdd', 'cwd', 'r95ptot', 'r99ptot', 'prcptot')

  indices.that.should.work <- list(tmax.indices, tmin.indices, c(tmax.indices, tmin.indices, tmax.tmin.indices), prec.indices, c(prec.indices, tmax.indices), c(prec.indices, tmin.indices), c(prec.indices, tmax.indices, tmin.indices, tmax.tmin.indices))
  
  for(i in 1:7) {
    include.tmax <- i %% 2
    include.tmin <- floor(i / 2) %% 2
    include.prec <- floor(i / 4) %% 2
    print(i)
    ci <- climdexInput.raw(if(include.tmax) tmax.dat else NULL,
                           if(include.tmin) tmin.dat else NULL,
                           if(include.prec) prec.dat else NULL,
                           if(include.tmax) tmax.dates else NULL,
                           if(include.tmin) tmin.dates else NULL,
                           if(include.prec) prec.dates else NULL,
                           base.range=c(1971, 2000))

    outbase.thresholds <- get.outofbase.quantiles(if(include.tmax) tmax.dat else NULL,
                           if(include.tmin) tmin.dat else NULL,
                           if(include.prec) prec.dat else NULL,
                           if(include.tmax) tmax.dates else NULL,
                           if(include.tmin) tmin.dates else NULL,
                           if(include.prec) prec.dates else NULL,
                           base.range=c(1971, 2000))

    ci.csv <- climdexInput.csv(if(include.tmax) "1018935_MAX_TEMP.csv" else NULL,
                               if(include.tmin) "1018935_MIN_TEMP.csv" else NULL,
                               if(include.prec) "1018935_ONE_DAY_PRECIPITATION.csv" else NULL,
                               data.columns=list(tmax="MAX_TEMP", tmin="MIN_TEMP", prec="ONE_DAY_PRECIPITATION"),
                               base.range=c(1971, 2000))

    indices.to.check.equals <- indices.that.should.work[[i]]
    indices.to.check.error <- all.indices[!(all.indices %in% indices.to.check.equals)]
    
    print("Checking success...")
    for(index in indices.to.check.equals) {
      print(index)
      fun <- match.fun(paste('climdex', index, sep="."))
      valid.result <- get(paste('ec.1018935', index, sep="."))
      checkEquals(valid.result, fun(ci), paste(index, "didn't match expected result: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
      checkEquals(valid.result, fun(ci.csv), paste(index, "didn't match expected result with CSV input: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
    }

    print("Checking error...")
    for(index in indices.to.check.error) {
      print(index)
      fun <- match.fun(paste('climdex', index, sep="."))
      checkException(fun(ci), paste(index, "didn't produce error when it should have: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
      checkException(fun(ci.csv), paste(index, "didn't produce error when it should have with CSV input: include.tmax is", include.tmax, ", include.tmin is ", include.tmin, ", include.prec is ", include.prec))
    }
  }
}

climdex.pcic.test.climdex.gsl <- function() {
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year", "jday")]), format="%Y %j", cal="gregorian")

  ## Load the data in.
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION, tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000), northern.hemisphere=FALSE)

  checkEquals(get('ec.1018935.gsl.sh'), climdex.gsl(ci))
}
