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
