library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"
source("../../climdex.pcic/R/climdex.r")

ClimVars.tmax <- file.path(wd, "1018935_MAX_TEMP.csv")
ClimVars.tmin <- file.path(wd, "1018935_MIN_TEMP.csv")
ClimVars.prec <- file.path(wd, "1018935_ONE_DAY_PRECIPITATION.csv")
datacol <- list(tmax = "MAX_TEMP", tmin = "MIN_TEMP", prec = "ONE_DAY_PRECIPITATION")

# Create climdex input objects for both northern and southern hemispheres.
ci.csv <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990))
ci.csv.sh <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990), northern.hemisphere = F)

# Returns the specific value of a climate variable at a specific index. Used in n or x, and Rxnday tests.
expected.result <- function(ci.csv, data, date.index, date.factors, freq, na.mask) {
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
# Generic function to get the expected results for N or X indices.
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
  date.index <- tapply(data, ci.csv@date.factors[[match.arg(freq)]], function(x) fun(x))
  date.factors <- unique(ci.csv@date.factors[[match.arg(freq)]])
  return(expected.result(ci.csv, data, date.index, date.factors, freq, na.mask))
}

# Test for TXx, TNn, TNx and TXn indices
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

# Generic function to get results for rx1day and rx5day indices.
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
  return(expected.result(ci.csv, data, date.index, date.factors, freq, na.mask))
}

# Test exact dates returned for Rx1day and Rx5day indices.
climdex.pcic.test.exact.date.rxnd.indices <- function() {
  test.indices <- c("rx1day", "rx5day")
  date.factors <- c("annual", "monthly", "seasonal")

  for (idx in test.indices) {
    for (freq in date.factors) {
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

# Find the longest consecutive true values in a boolean array.
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

# Function to get expected spell boundaries for CDD and CWD indices.
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
# Generic to compare the expected and climdex-calculated results for the spell tests.
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

# Test cdd and cwd spells with example data set.
climdex.pcic.test.spell.boundaries <- function() {
  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    if (idx == "cdd") {
      result <- climdex.cdd(ci.csv, spells.can.span.years = F, as.df = TRUE)
    } else {
      result <- climdex.cwd(ci.csv, spells.can.span.years = F, as.df = TRUE)
    }

    expected <- get.spell.bounds(ci.csv, idx)
    check.spell.results(expected, result, idx)
  }
}

# Test for spell that spans multiple years.
climdex.pcic.test.multi.year.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")

  test.prec.data.dry <- rep(0, length(test.dates))
  ci.dry <- climdexInput.raw(prec = test.prec.data.dry, prec.dates = test.dates)

  test.prec.data.wet <- rep(2, length(test.dates))
  ci.wet <- climdexInput.raw(prec = test.prec.data.wet, prec.dates = test.dates)

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
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

# Test for zero-length spell duration.
climdex.pcic.test.no.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.prec.data.dry <- rep(0, length(test.dates))
  test.prec.data.wet <- rep(1, length(test.dates))

  ci.wet <- climdexInput.raw(prec = test.prec.data.wet, prec.dates = test.dates)
  ci.dry <- climdexInput.raw(prec = test.prec.data.dry, prec.dates = test.dates)

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
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

# Test function for spell boundaries with NA masks.
climdex.pcic.test.na.masks.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.dates.factor <- factor(format(test.dates, "%Y"))
  test.prec.data.ran <- c(sample(0:2, length(test.dates), replace = TRUE))
  ci.ran <- climdexInput.raw(prec = test.prec.data.ran, prec.dates = test.dates)
  ci.ran@namasks$annual$prec <- rep(c(NA, 1), length.out = length(unique(test.dates.factor)))

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    if (idx == "cdd") {
      result <- climdex.cdd(ci.ran, spells.can.span.years = F, as.df = TRUE)
    } else {
      result <- climdex.cwd(ci.ran, spells.can.span.years = F, as.df = TRUE)
    }
    
    expected <- get.spell.bounds(ci.ran, idx)
    check.spell.results(expected, result, idx)
  }
}

# Return start or end of the GSL index as PCICt object in %Y-%m-%d format
gsl.test.ymd <- function(year, cal, doy, northern.hemisphere) {
  origin <- ifelse(northern.hemisphere, paste(year, "01", "01", sep = "-"), paste(year, "07", "01", sep = "-"))
  origin.pcict <- as.PCICt(origin, cal)
  seconds.per.day <- 86400
  doy.pcict <- origin.pcict + (doy) * seconds.per.day
  ymd <- as.PCICt(doy.pcict, cal = cal, format = "%Y-%m-%d")
  return(format(ymd, "%Y-%m-%d"))
}

# Helper for finding the GSL bounds.
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

# Calculate the Expected GSL bounds.
expected.gsl <- function(ci, as.df) {
  cal <- attr(ci@dates, "cal")
  n.h <- ci@northern.hemisphere

  gsl.output <- lapply(unique(ci@date.factors$annual), function(year) {
    year.length <- length(ci@dates[ci@date.factors$annual == year])
    year.idx <- ci@date.factors$annual == year
    next.year <- as.numeric(as.character(year)) + 1
    leap.year <- year.length > 365
    next.year.is.leap <- length(ci@dates[ci@date.factors$annual == next.year]) > 365

    if (!n.h) {
      sy <- paste(year, "07", "01", sep = "-")
      sy.pcict <- as.PCICt(sy, cal)
      seconds.per.day <- 86400
      ey.pcict <- sy.pcict + ifelse(next.year.is.leap, year.length + 1, ifelse(leap.year, year.length - 1, year.length)) * seconds.per.day
      date_indices <- which(ci@dates >= sy.pcict & ci@dates < ey.pcict)

      tavg <- ci@data$tavg[date_indices]
    } else {
      tavg <- ci@data$tavg[year.idx]
    }

    midpoint <- ceiling(length(tavg) / 2)
    if (next.year.is.leap) midpoint <- midpoint + 1

    first.half <- tavg[1:(midpoint - 1)]
    second.half <- tavg[midpoint:length(tavg)]

    start.idx <- find.repeated(first.half, 6, ">")
    end.idx <- find.repeated(second.half, 6, "<") - 1

    start.result <- ifelse(is.na(start.idx), NA_character_, gsl.test.ymd(year, cal, start.idx - 1, n.h))

    if (is.na(start.idx)) {
      end.result <- NA_character_
      duration <- 0
    } else if (is.na(end.idx)) {
      end.result <- gsl.test.ymd(year, cal, length(tavg), n.h)
      duration <- length(tavg) + 1 - start.idx
    } else {
      end.result <- gsl.test.ymd(year, cal, (ifelse(end.idx > 1, midpoint + end.idx - 1, ifelse(!n.h && next.year.is.leap, midpoint, midpoint + 1))), n.h)
      duration <- (midpoint + end.idx) - start.idx
    }

    list(start = start.result, end = end.result, duration = duration)
  })
  starts <- lapply(gsl.output, `[[`, "start")
  ends <- lapply(gsl.output, `[[`, "end")
  season <- lapply(gsl.output, `[[`, "duration")

  ex.r <- data.frame(start = unlist(starts), end = unlist(ends), sl = unlist(season))
  ex.r$sl <- ex.r$sl * ci@namasks$annual$tavg
  ex.r$start[is.na(ex.r$sl)] <- NA
  ex.r$end[is.na(ex.r$sl)] <- NA

  if (as.df) {
    ex.r <- data.frame(start = unlist(starts), end = unlist(ends), sl = unlist(season))

    if (!n.h) {
      ci@namasks$annual$tavg[length(ci@namasks$annual$tavg)] <- NA
      dates.POSIXlt <- as.POSIXlt(ci@dates)
      years <- dates.POSIXlt$year + 1900
      months <- dates.POSIXlt$mon + 1

      valid.years <- range(years)
      years.gsl <- years - floor((12 - months) / 6)

      inset <- years.gsl >= valid.years[1]
      gsl.factor <- factor(years.gsl[inset])
      gsl.factor.monthly <- factor(paste(years.gsl[inset], months[inset], sep = "-"))
      gsl.yearmonth.factor <- unlist(strsplit(levels(gsl.factor.monthly), "-"))[(0:(nlevels(gsl.factor.monthly) - 1)) * 2 + 1]
      gsl.temp.data <- ci@data$tavg[inset]
      namask.gsl.monthly <- get.na.mask(gsl.temp.data, gsl.factor.monthly, ci@max.missing.days["annual"])
      ci@namasks$annual$tavg <- get.na.mask(gsl.temp.data, gsl.factor, ci@max.missing.days["annual"]) * as.numeric(tapply(namask.gsl.monthly, gsl.yearmonth.factor, prod))
    }

    ex.r$sl <- ex.r$sl * ci@namasks$annual$tavg
    ex.r$start[is.na(ex.r$sl)] <- NA
    ex.r$end[is.na(ex.r$sl)] <- NA
  }

  return(ex.r)
}
# Generic GSL test function comparing expected with climdex-calculated results.
test_gsl <- function(ci, test_name) {
  expected <- expected.gsl(ci, as.df = TRUE)
  result <- climdex.gsl(ci, "GSL", as.df = TRUE)

  for (i in seq_along(result$start)) {
    if (test_name == "climdex.pcic.test.na.masks.gsl") {
      # Only compare NA rows
      if (i %% 2 == 0) {
        next
      }
    }
    checkIdentical(expected$start[i], result$start[i], paste("Start of GSL does not match Expected:", as.character(expected$start[i]), " Result: ", result$start[i], " (", test_name, ")"))
    checkIdentical(expected$end[i], result$end[i], paste("End of GSL does not match Expected:", as.character(expected$end[i]), " Result: ", result$end[i], " (", test_name, ")"))
  }
}

# Test GSL calculation with example data.
climdex.pcic.test.gsl <- function() {
  test_gsl(ci.csv, "climdex.pcic.test.gsl")
}

# Test for a GSL that spans multiple years.
climdex.pcic.test.multi.year.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.tavg.data <- c(rep(4, cal - 15), rep(6, length(test.dates) - (cal - 15)))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  test_gsl(ci.gsl, "climdex.pcic.test.multi.year.gsl")
}

# Test for GSL that never starts.
climdex.pcic.test.no.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1964-12-31", cal = cal), by = "days")
  test.tavg.data <- c(rep(4, length(test.dates)))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  test_gsl(ci.gsl, "climdex.pcic.test.no.gsl")
}

# Test that the NA masks are applied to results.
climdex.pcic.test.na.masks.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1969-12-31", cal = cal), by = "days")
  test.dates.factor <- factor(format(test.dates, "%Y"))
  test.tavg.data <- c(sample(1:10, length(test.dates), replace = TRUE))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  ci.gsl@namasks$annual$tavg <- rep(c(NA, 1), length.out = length(unique(test.dates.factor)))
  test_gsl(ci.gsl, "climdex.pcic.test.na.masks.gsl")
}

# Test GSL in the southern hemisphere.
climdex.pcic.test.gsl.southern.hemisphere <- function() {
  test_gsl(ci.csv.sh, "climdex.pcic.test.gsl.southern.hemisphere")
}
