library(RUnit)

source('../climdex.r', chdir=T)

make.temperature.fixture <- function() {
  list(
       list(temp=-5:5,         as.factor(rep('2010', 11))),
       list(temp=c(-5:5,5:-5), as.factor(c(rep('2010', 11), rep('2011', 11)))),
       ## All are negative
       list(temp=rep(-10, 20), as.factor(c(rep('2010', 10), rep('2011', 10)))),
       ## All are positive
       list(temp=rep(10, 20), as.factor(c(rep('2010', 10), rep('2011', 10))))
       )
}

mk.rv <- function(a, names)
  array(a, dimnames=list(names))

test.can.load.data <- function() {
  station <- '1098D90'
  data.dir <- '/home/data/projects/data_cleanup/CDCD_2007/new_data/'
  vars <- c(tmax='MAX_TEMP', tmin='MIN_TEMP', prec='ONE_DAY_PRECIPITATION')
  data.files <- file.path(data.dir, vars, paste(station, '_', vars, '.csv', sep=''))
  args <- append(data.files, list(data.column=as.list(vars)))
  clim.in <- do.call(climdexInput, args)
  checkTrue(inherits(clim.in, 'climdexInput'))
}

test.CSDI <- function() {
  f <- CSDI
  cases <- list(list(args=list(n=3, v=c(F, T, T, T, T, F), years=as.factor(rep('2008', 6))),
                     expected=mk.rv(4, '2008')),
                list(args=list(n=3, v=c(F, T, T, T, T, F), years=as.factor(c(rep('2008', 2), rep('2009', 4)))),
                     expected=mk.rv(c(1, 3), c('2008', '2009'))),
                list(args=list(n=3, v=c(F, rep(T, 4), rep(F, 2), rep(T, 4)), years=as.factor(c(rep('2008', 6), rep('2009', 5)))),
                     expected=mk.rv(c(4, 4), c('2008', '2009'))),
                list(args=list(n=0, v=c(T, F, T, F, T, F), years=as.factor(rep('2008', 6))),
                     expected=mk.rv(3, '2008'))
                )
  for (case in cases) {
    checkEquals(do.call(f, case$args), case$expected)
  }
  checkException(f(n=-1, v=rep(T, 10)))
}

test.number.days.below.threshold <- function() {
  f <- number.days.below.threshold
  fix <- make.temperature.fixture()
  thresh <- 0
  expected <- list(mk.rv(5, '2010'),
                   mk.rv(c(5, 5), c('2010', '2011')),
                   mk.rv(c(10, 10), c('2010', '2011')),
                   mk.rv(c(0, 0), c('2010', '2011'))
                   )
  do.check <- function(args, y) {
    args <- append(append(args, f, after=0), thresh)
    cl <- as.call(args)
    checkEquals(eval(cl), y)
  }
  mapply(do.check, fix, expected)

  ## Should raise exception
  ex.checks <- list(list(f, NULL),
                    ## Lengths don't match
                    list(f, 1:5, as.factor(rep('2010', 4)), 0),
                    ## Not numeric
                    list(f, 'blah stuff', as.factor('2010'), 0)
                    )
  rv <- lapply(ex.checks, as.call)
  lapply(rv, checkException)
}

test.number.days.over.threshold <- function() {
  f <- number.days.over.threshold
  fix <- make.temperature.fixture()
  thresh <- 0
  expected <- list(mk.rv(5, '2010'),
                   mk.rv(c(5, 5), c('2010', '2011')),
                   mk.rv(c(0, 0), c('2010', '2011')),
                   mk.rv(c(10, 10), c('2010', '2011'))
                   )
  do.check <- function(args, y) {
    args <- append(append(args, f, after=0), thresh)
    cl <- as.call(args)
    checkEquals(eval(cl), y)
  }
  mapply(do.check, fix, expected)

  ## Should raise exception
  ex.checks <- list(list(f, NULL),
                    ## Lengths don't match
                    list(f, 1:5, as.factor(rep('2010', 4)), 0),
                    ## Not numeric
                    list(f, 'blah stuff', as.factor('2010'), 0)
                    )
  rv <- lapply(ex.checks, as.call)
  lapply(rv, checkException)
}

## Takes
check.one.case <- function(case, f) {
  args <- append(case[- which(names(case) == 'expected')], f, after=0)
  cl <- as.call(args)
  print(eval(cl))
  checkEquals(eval(cl), case$expected)
}

test.growing.season.length <- function() {
  f <- growing.season.length

  fac <- as.factor(c(rep('2010', 366), rep('2011', 365), rep('2012', 365)))

  cases <- list()

  ## First and easiest case: a 6 day season beginning well into the second half of the year
  twenty.o.nine <- seq(as.POSIXct("2009/01/01", tz="GMT"), by="day", length.out=365)
  x <- my.ts(rep(0, 365), twenty.o.nine)
  x[seq(as.POSIXct("2009/08/01", tz="GMT"), by="day", length.out=6)] <- rep(5.1, 6)
  expected <- mk.rv(6, '2009')
  cases <- append(cases, list(list(x, as.factor(rep('2009', 365)), expected=expected)))

  ## Season starts at the beginning of the year, ends after July 1
  x <- my.ts(rep(10, 365), twenty.o.nine)
  fac <- as.factor(rep('2009', 365))
  x[seq(as.POSIXct("2009/07/02", tz="GMT"), by="day", length.out=6)] <- rep(0)
  expected <- mk.rv(as.POSIXlt("2009/07/01", tz="GMT")$yday, '2009')
  cases <- append(cases, list(list(x, as.factor(rep('2009', 365)), expected=expected)))

  ## Simple case: 1,2,3 day seasons starting right after July 1
  ## wedged in between a 6 day series of 5.1 and zeros
  x <- my.ts(rep(0, 365*3+1), seq(as.POSIXct("2010/01/01", tz="GMT"), by="day", length.out=365*3+1))
  for (year in 2010:2012) {
    x[seq(as.POSIXct(paste(year, "/07/01", sep=""), tz="GMT"), length.out=6, by=as.difftime(-1, units="days"))] <- rep(5.1, 6)
    x[seq(as.POSIXct(paste(year, "/07/02", sep=""), tz="GMT"), length.out=year-2009, by=as.difftime(1, units="days"))] <- rep(10, year-2009)
  }
  expected <- mk.rv(7:9, c('2010', '2011', '2012'))
  cases <- append(cases, list(list(x, fac, expected=expected)))

  lapply(cases, check.one.case, f)
}

## These are simple enough definitions, that I don't think they need testing
test.max.min.daily.temp <- function() {
}



## Utility timeseries class
## Carries around a POSIXct vector as an attribute that describes the
## timeseries.  Allows subsetting and subset replacement by passing in
## POSIX types as indicies
my.ts <- function(x, t) {
  stopifnot(length(x) == length(t))
  stopifnot(inherits(t, "POSIXt"))
  class(x) <- append(after=0, class(x), "my.ts")
  attr(x, 'time') <- t
  x
}

'[.my.ts' <- function(obj, ti) {
  if (is.numeric(ti)) {
    f <- '['
    x <- eval(call(f, as.numeric(obj), ti))
    return(my.ts(x, attr(obj, 'time')[ti]))
  } else if (inherits(ti, 'POSIXt')) {
    ## This will be really slow for long arrays... but it's simple
    i <- sapply(ti, function(x) {which(attr(obj, 'time') == x)})
    return(obj[i])
  }

}

'[<-.my.ts' <- function(obj, ti, value) {
  if (is.numeric(ti)) {
    f <- '[<-'
    x <- eval(call(f, as.numeric(obj), ti, value))
    return(my.ts(x, attr(obj, 'time')))
  } else if (inherits(ti, 'POSIXt')) {
    ## This will be really slow for long arrays... but it's simple
    i <- sapply(ti, function(x) {which(attr(obj, 'time') == x)})
    obj[i] <- value
    return(obj)
  }
}
