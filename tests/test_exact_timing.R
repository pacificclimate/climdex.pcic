library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"
source("../../climdex.pcic/R/climdex.r")

ClimVars.tmax <- file.path(wd, "1018935_MAX_TEMP.csv")
ClimVars.tmin <- file.path(wd, "1018935_MIN_TEMP.csv")
ClimVars.prec <- file.path(wd, "1018935_ONE_DAY_PRECIPITATION.csv")


datacol <- list(tmax = "MAX_TEMP", tmin = "MIN_TEMP", prec = "ONE_DAY_PRECIPITATION")
ci.csv <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990))
ci.csv.sh <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990), northern.hemisphere = F)

expected.result <- function(ci.csv, data, date.index, date.factors, freq, fun, offset, na.mask) {
  ex.r <- c()
  ex.r <- lapply(seq_along(date.factors), function(i) {
    factor.name <- names(date.index[i])
    dof <- date.index[[i]]
    dates <- ci.csv@dates[ci.csv@date.factors[[freq]] == factor.name]
    if (!is.na(na.mask[i])) {
      dates[dof]
    } else {
      NA
    }
  })
  return((ex.r))
}


get.n.or.x.result <- function(idx, ci.csv, freq = c("monthly", "annual", "seasonal")) {
  data <- if (substr(idx, 2, 2) == "x") {
    ci.csv@data$tmax
  } else {
    ci.csv@data$tmin
  }
  na.mask <- if (substr(idx, 2, 2) == "x") {
    ci.csv@namasks[[freq]]$tmax
  } else {
    ci.csv@namasks[[freq]]$tmin
  }
  fun <- ifelse(substr(idx, 3, 3) == "x", which.max, which.min)
  offset <- 0
  date.index <- tapply(data, ci.csv@date.factors[[match.arg(freq)]], function(x) fun(x))
  date.factors <- unique(ci.csv@date.factors[[match.arg(freq)]])
  return(expected.result(ci.csv, data, date.index, date.factors, freq, fun, offset, na.mask))
}

# Exact date:
climdex.pcic.test.exact.date.n.or.x.indices <- function() {
  test.indices <- c("txx", "tnn", "tnx", "txn")
  date.factors <- c("annual", "monthly", "seasonal")

  for (idx in test.indices) {
    data <- if (substr(idx, 2, 2) == "x") {
      ci.csv@data$tmax
    } else {
      ci.csv@data$tmin
    }

    for (freq in date.factors) {
      print(paste(idx, freq))
      fun <- paste("climdex", idx, sep = ".")
      result <- do.call(fun, list(ci.csv, freq = freq, as.df = TRUE))
      expected <- get.n.or.x.result(idx, ci.csv, freq)
      result$ymd <- as.character(result$ymd)
      for (i in seq_along(expected)) {
        expected.val <- data[ci.csv@dates %in% as.character(expected[[i]])]
        expected.val <- ifelse(length(expected.val) == 0, NA, expected.val)
        checkIdentical(as.character(expected[[i]]), result$ymd[i])
        checkEqualsNumeric(as.numeric(expected.val), as.numeric(result$val[i]))
      }
    }
  }
}


get.Rxnday.result <- function(idx, ci.csv, freq = c("monthly", "annual", "seasonal"), ndays, center.mean.on.last.day) {
  data <- ci.csv@data$prec
  na.mask <- ci.csv@namasks[[freq]]$prec
  fun <- which.max
  if (as.numeric(substr(idx, 3, 3)) == 5) {
    data[is.na(data)] <- 0
    prec.runsum <- running.mean(data, ndays)
    prec.runsum[is.na(prec.runsum)] <- 0

    if (center.mean.on.last.day) {
      k2 <- ndays %/% 2
      prec.runsum <- c(rep(0, k2), prec.runsum[1:(length(prec.runsum) - k2)])
    }
    data <- prec.runsum
  }
  date.index <- tapply(data, ci.csv@date.factors[[match.arg(freq)]], function(x) fun(x))
  date.factors <- unique(ci.csv@date.factors[[match.arg(freq)]])
  return(expected.result(ci.csv, data, date.index, date.factors, freq, fun, offset, na.mask))
}


climdex.pcic.test.exact.date.rxnd.indices <- function() {
  test.indices <- c("rx1day", "rx5day")
  date.factors <- c("annual", "monthly", "seasonal")

  for (idx in test.indices) {
    for (freq in date.factors) {
      print(paste(idx, freq))
      fun <- paste("climdex", idx, sep = ".")
      result <- do.call(fun, list(ci.csv, freq = freq, as.df = TRUE))
      ndays <- as.numeric(substr(idx, 3, 3))
      center.mean.on.last.day <- FALSE
      expected <- get.Rxnday.result(idx, ci.csv, freq, ndays, center.mean.on.last.day)
      result$ymd <- as.character(result$ymd)

      for (i in seq_along(result$ymd)) {
        if (!is.na(result$ymd[i])) {
          if (ndays == 5) {
            window.start <- expected[[i]] - 2 * 86400
            window.end <- expected[[i]] + 2 * 86400
            expected.value <- sum(ci.csv@data$prec[ci.csv@dates >= window.start & ci.csv@dates <= window.end], na.rm = TRUE)
          } else {
            expected.value <- ci.csv@data$prec[ci.csv@dates == expected[[i]]]
          }
        } else {
          expected.value <- NA
          checkTrue(is.na(expected[[i]]) && is.na(result$val[i]))
        }
        checkEqualsNumeric(as.numeric(expected.value), as.numeric(result$val[i]), tolerance = 0.01)
        checkIdentical(as.character(expected[[i]]), result$ymd[i])
      }
    }
  }
}


find.longest.consecutive.true <- function(bool.array) {
  max.length <- current.length <- start.index <- max.start.index <- 0

  for (i in seq_along(bool.array)) {
    if (!is.na(bool.array[i]) && bool.array[i]) {
      current.length <- current.length + 1
      if (current.length == 1) start.index <- i
      if (start.index) {
        if (current.length > max.length) {
          max.length <- current.length
          max.start.index <- start.index
        }
      }
    } else {
      start.index <- current.length <- 0
    }
  }

  data.frame(start = max.start.index, end = max.start.index + max.length - 1, duration = max.length)
}


get.spell.bounds <- function(ci, idx) {
  spell.boundary <- lapply(unique(ci@date.factors$annual), function(year) {
    if (idx == "cdd") {
      bool.array <- ci@data$prec[ci@date.factors$annual == year] < 1
    } else {
      bool.array <- ci@data$prec[ci@date.factors$annual == year] >= 1
    }
    spell <- find.longest.consecutive.true(bool.array)
    dates <- ci@dates[ci@date.factors$annual == year]
    if (spell$duration > 0) {
      spell.bounds <- data.frame(start = dates[spell$start], end = dates[spell$end], duration = ifelse(is.na(ci@namasks$annual$prec[year]), NA, spell$duration))
    } else {
      spell.bounds <- data.frame(start = NA, end = NA, duration = 0)
    }

    return(spell.bounds)
  })
  return(spell.boundary)
}

check.spell.results <- function(expected, result, idx) {
  for (i in seq_along(result$start)) {
    if (is.na(expected[[i]]$duration)) {
      expected[[i]]$start <- NA
      expected[[i]]$end <- NA
    }
    checkIdentical(as.character(expected[[i]]$start), result$start[i], paste("Start of spells for index: ", idx, " do not agree: Expected:", as.character(expected[[i]]$start), " Result: ", result$start[i]))
    checkIdentical(as.character(expected[[i]]$end), result$end[i])
  }
}


climdex.pcic.test.spell.boundaries <- function() {
  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    print(paste(idx))

    if (idx == "cdd") {
      result <- climdex.cdd(ci.csv, spells.can.span.years = F, as.df = TRUE)
    } else {
      result <- climdex.cwd(ci.csv, spells.can.span.years = F, as.df = TRUE)
    }

    expected <- get.spell.bounds(ci.csv, idx)
    check.spell.results(expected, result, idx)
  }
}


climdex.pcic.test.multi.year.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")

  test.prec.data.dry <- rep(0, length(test.dates))
  ci.dry <- climdexInput.raw(prec = test.prec.data.dry, prec.dates = test.dates)

  test.prec.data.wet <- rep(2, length(test.dates))
  ci.wet <- climdexInput.raw(prec = test.prec.data.wet, prec.dates = test.dates)

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    print(paste(idx))

    if (idx == "cdd") {
      result <- climdex.cdd(ci.dry, spells.can.span.years = F, as.df = TRUE)
      ci <- ci.dry
    } else {
      result <- climdex.cwd(ci.wet, spells.can.span.years = F, as.df = TRUE)
      ci <- ci.wet
    }

    expected <- get.spell.bounds(ci, idx)
    check.spell.results(expected, result, idx)
  }
}


climdex.pcic.test.no.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.prec.data.dry <- rep(0, length(test.dates))
  test.prec.data.wet <- rep(1, length(test.dates))

  ci.wet <- climdexInput.raw(prec = test.prec.data.wet, prec.dates = test.dates)
  ci.dry <- climdexInput.raw(prec = test.prec.data.dry, prec.dates = test.dates)

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    print(paste(idx))

    if (idx == "cdd") {
      result <- climdex.cdd(ci.wet, spells.can.span.years = F, as.df = TRUE)
      ci <- ci.wet
    } else {
      result <- climdex.cwd(ci.dry, spells.can.span.years = F, as.df = TRUE)
      ci <- ci.dry
    }

    expected <- get.spell.bounds(ci, idx)
    check.spell.results(expected, result, idx)
  }
}


climdex.pcic.test.na.masks.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.dates.factor <- factor(format(test.dates, "%Y"))
  test.prec.data.ran <- c(sample(0:2, length(test.dates), replace = TRUE))
  ci.ran <- climdexInput.raw(prec = test.prec.data.ran, prec.dates = test.dates)
  ci.ran@namasks$annual$prec <- rep(c(NA, 1), length.out = length(unique(test.dates.factor)))

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    print(paste(idx))

    if (idx == "cdd") {
      result <- climdex.cdd(ci.ran, spells.can.span.years = F, as.df = TRUE)
    } else {
      result <- climdex.cwd(ci.ran, spells.can.span.years = F, as.df = TRUE)
    }
    expected <- get.spell.bounds(ci.ran, idx)
    check.spell.results(expected, result, idx)
  }
}


gsl.test.ymd <- function(year, cal, doy, northern.hemisphere) {
  origin <- ifelse(northern.hemisphere, paste(year, "01", "01", sep = "-"), paste(year, "07", "01", sep = "-"))
  origin.pcict <- as.PCICt(origin, cal)
  seconds.per.day <- 86400
  doy.pcict <- origin.pcict + (doy - 1) * seconds.per.day
  ymd <- as.PCICt(doy.pcict, cal = cal, format = "%Y-%m-%d")
  return(format(ymd, "%Y-%m-%d"))
}


find.repeated <- function(tavg, repetition, op) {
  count <- 0
  index <- 0
  set.flag <- F
  for (i in seq_along(tavg)) {
    if (!is.na(tavg[i])) {
      if ((op == ">" && tavg[i] > 5) || (op == "<" && tavg[i] < 5)) {
        count <- count + 1
      } else {
        count <- 0
      }

      if (count >= repetition) {
        index <- i - repetition + 1
        set.flag <- T
        break
      }
    } else {
      count <- 0
    }
  }
  if (set.flag) {
    return(index)
  } else {
    return(NA)
  }
}


expected.gsl <- function(ci, as.df) {
  cal <- attr(ci@dates, "cal")
  n.h <- ci@northern.hemisphere
  gsl.output <- lapply(unique(ci@date.factors$annual), function(year) {
    year.length <- length(ci@dates[ci@date.factors$annual == year])
    year.idx <- ci@date.factors$annual == year

    if (!n.h) {
      sy <- paste(year, "07", "01", sep = "-")
      sy.pcict <- as.PCICt(sy, cal)
      seconds.per.day <- 86400
      ey.pcict <- sy.pcict + year.length * seconds.per.day
      date_indices <- which(ci@dates >= sy.pcict & ci@dates < ey.pcict)
      tavg <- ci@data$tavg[date_indices]
    } else {
      tavg <- ci@data$tavg[year.idx]
    }

    midpoint <- length(tavg) %/% 2
    first.half <- tavg[1:midpoint - 1]
    second.half <- tavg[(midpoint):length(tavg)]
    op.gt <- ">"
    op.lt <- "<"

    start.idx <- find.repeated(first.half, 6, op.gt)
    end.idx <- find.repeated(second.half, 6, op.lt) - 1
    start.result <- gsl.test.ymd(year, cal, start.idx, n.h)

    if (is.na(start.idx)) {
      start.result <- NA_character_
      end.result <- NA_character_
      duration <- 0
    } else if (is.na(end.idx)) {
      end.result <- gsl.test.ymd(year, cal, length(tavg) + 1, n.h)
      duration <- length(tavg) + 1 - start.idx
    } else {
      end.result <- gsl.test.ymd(year, cal, midpoint + end.idx, n.h)
      duration <- midpoint + end.idx - start.idx
    }
    list(start = start.result, end = end.result, duration = duration)
  })

  starts <- lapply(gsl.output, `[[`, "start")
  ends <- lapply(gsl.output, `[[`, "end")
  season <- lapply(gsl.output, `[[`, "duration")

  if (as.df) {
    ex.r <- data.frame(start = unlist(starts), end = unlist(ends), sl = unlist(season))
    if (!n.h) ci@namasks$annual$tavg[length(ci@namasks$annual$tavg)] <- NA
    ex.r$sl <- ex.r$sl * ci@namasks$annual$tavg
    ex.r$start[is.na(ex.r$sl)] <- NA
    ex.r$end[is.na(ex.r$sl)] <- NA
  }
  return(ex.r)
}

# Generic GSL test function
test_gsl <- function(ci, test_name) {
  expected <- expected.gsl(ci, as.df = TRUE)
  result <- climdex.gsl(ci, "GSL", as.df = TRUE)

  for (i in seq_along(result$start)) {
    checkIdentical(expected$start[i], result$start[i], paste("Start of GSL does not match Expected:", as.character(expected$start[i]), " Result: ", result$start[i], " (", test_name, ")"))
    checkIdentical(expected$end[i], result$end[i], paste("End of GSL does not match Expected:", as.character(expected$end[i]), " Result: ", result$end[i], " (", test_name, ")"))
  }
}

climdex.pcic.test.gsl <- function() {
  test_gsl(ci.csv, "climdex.pcic.test.gsl")
}


climdex.pcic.test.multi.year.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.tavg.data <- c(rep(4, cal - 15), rep(6, length(test.dates) - (cal - 15)))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  test_gsl(ci.gsl, "climdex.pcic.test.multi.year.gsl")
}


climdex.pcic.test.no.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1964-12-31", cal = cal), by = "days")
  test.tavg.data <- c(rep(4, length(test.dates)))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  test_gsl(ci.gsl, "climdex.pcic.test.no.gsl")
}


climdex.pcic.test.na.masks.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1969-12-31", cal = cal), by = "days")
  test.dates.factor <- factor(format(test.dates, "%Y"))
  test.tavg.data <- c(sample(1:10, length(test.dates), replace = TRUE))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  ci.gsl@namasks$annual$tavg <- rep(c(NA, 1), length.out = length(unique(test.dates.factor)))
  test_gsl(ci.gsl, "climdex.pcic.test.na.masks.gsl")
}
