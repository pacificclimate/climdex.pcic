library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"

ClimVars.tmax <- file.path(wd, "1018935_MAX_TEMP.csv")
ClimVars.tmin <- file.path(wd, "1018935_MIN_TEMP.csv")
ClimVars.prec <- file.path(wd, "1018935_ONE_DAY_PRECIPITATION.csv")
datacol <- list(tmax = "MAX_TEMP", tmin = "MIN_TEMP", prec = "ONE_DAY_PRECIPITATION")

# Create climdex input objects for both northern and southern hemispheres.
ci.csv <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990))
ci.csv.sh <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990), northern.hemisphere = F)

expected.exact.date <- function(ci.csv, data, factor.extremes, date.factors, freq, na.mask) {
  ex.r <- c()
  ex.r <- lapply(seq_along(date.factors), function(i) {
    factor.name <- names(factor.extremes[i])
    idx.of.extreme <- factor.extremes[[i]]
    dates.in.factor <- ci.csv@dates[ci.csv@date.factors[[freq]] == factor.name]

    if (!is.na(na.mask[i])) {
      dates.in.factor[idx.of.extreme]
    } else {
      NA_character_
    }
  })

  return(ex.r)
}

get.data.for.idx <- function(ci, idx) {
  var <- climdex.min.max.idx.list[[idx]]$var
  ci@data[[var]]
}

are.not.all.na <- function(x,r) {
  checkTrue(any(!is.na(x)))
  checkTrue(any(!is.na(r)))
}

get.data.for.idx <- function(ci, idx) {
  if (substr(idx, 2, 2) == "x") {
    ci@data$tmax
  } else {
    ci@data$tmin
  }
}

are.not.all.na <- function(x,r) {
  checkTrue(any(!is.na(x)))
  checkTrue(any(!is.na(r)))
}
# Generic function to get the expected results for N or X indices.
get.n.or.x.result <- function(idx, ci.csv, freq = c("monthly", "annual", "seasonal")) {
  data <- get.data.for.idx(ci.csv, idx)
  na.mask <- if (substr(idx, 2, 2) == "x") {
    ci.csv@namasks[[freq]]$tmax
  } else {
    ci.csv@namasks[[freq]]$tmin
  }
  fun <- ifelse(substr(idx, 3, 3) == "x", which.max, which.min)
  factor.extremes <- tapply(data, ci.csv@date.factors[[match.arg(freq)]], fun)
  date.factors <- unique(ci.csv@date.factors[[match.arg(freq)]])
  return(expected.exact.date(ci.csv, data, factor.extremes, date.factors, freq, na.mask))
}

#Custom checkEqual to compare expected and climdex results.
is.almost.equal <- function(x, r, tolerance = 0.01) {
  if (is.na(x) && is.na(r)) {
    return(TRUE)
  } else if (is.na(x) || is.na(r)) {
    return(FALSE)
  } else {
    return(abs(x - r) < tolerance)
  }
}


# Test for TXx, TNn, TNx and TXn indices.
climdex.pcic.test.exact.date.n.or.x.indices <- function() {
  test.indices <- names(climdex.min.max.idx.list)[grepl("t", names(climdex.min.max.idx.list))]
  date.factors <- c("annual", "monthly", "seasonal")

  for (idx in test.indices) {
    data <- get.data.for.idx(ci.csv, idx)

    for (freq in date.factors) {
      fun <- paste("climdex", idx, sep = ".")
      result <- do.call(fun, list(ci.csv, freq = freq, include.exact.dates = TRUE))
      expected <- get.n.or.x.result(idx, ci.csv, freq)
      result$ymd <- as.character(result$ymd)
      checkIdentical(length(expected), nrow(result), paste("Lengths differ. Expected:", length(expected),"Result:", nrow(result)))
      are.not.all.na(expected, result$ymd)
       for (i in seq_along(expected)) {
        expected.val <- data[ci.csv@dates == expected[[i]]]
        expected.val <- ifelse(length(expected.val) == 0, NA, expected.val)
        checkIdentical(as.character(expected[[i]]), as.character(result$ymd[i]), 
                       paste("Idx:", idx, "Expected: ", as.character(expected[[i]]), "Result: ", as.character(result$ymd[i])))
        checkTrue(is.almost.equal(as.numeric(expected.val), as.numeric(result$val[i])), 
                  msg = paste("Idx:", idx, "Expected: ", as.numeric(expected.val), "Result: ", as.numeric(result$val[i])))
        
      }
    }
  }
}

# Boundary test that checks if TXx and TNn occur on the last day of the year,
# and TXn and TNx occur on the first day of the year.
climdex.pcic.test.n.or.x.dates.at.end.of.year <- function() {
  test.indices <- names(climdex.min.max.idx.list)[grepl("t", names(climdex.min.max.idx.list))]
  date.factors <- c("annual", "monthly", "seasonal")
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1961-12-31", cal = cal), by = "days")
  test.tmax <- c(rep(0, length(test.dates) - 1), 1)
  test.tmin <- c(rep(1, length(test.dates) - 1), 0)
  ci.nx.eoy <- climdexInput.raw(tmax = test.tmax, tmin = test.tmin, tmax.dates = test.dates, tmin.dates = test.dates)

  for (idx in test.indices) {
    data <- get.data.for.idx(ci.nx.eoy, idx)

    for (freq in date.factors) {
      fun <- paste("climdex", idx, sep = ".")

      result <- do.call(fun, list(ci.nx.eoy, freq = freq, include.exact.dates = TRUE))
      expected <- get.n.or.x.result(idx, ci.nx.eoy, freq)
      result$ymd <- as.character(result$ymd)
      checkIdentical(length(expected), nrow(result), paste("Lengths differ. Expected:", length(expected),"Result:", nrow(result)))
      are.not.all.na(expected, result$ymd)
      for (i in seq_along(expected)) {
        expected.val <- data[ci.nx.eoy@dates == expected[[i]]]
        expected.val <- ifelse(length(expected.val) == 0, NA, expected.val)
        checkIdentical(as.character(expected[[i]]), result$ymd[i], paste("Idx:", idx, "Expected: ", as.character(expected[[i]]), "Result: ", as.character(result$ymd[i])))
        checkTrue(is.almost.equal(as.numeric(expected.val), as.numeric(result$val[i])), 
                  msg = paste("Idx:", idx, "Expected: ", as.numeric(expected.val), "Result: ", as.numeric(result$val[i])))
      }
    }
  }
}

# Check extreme dates in winter season.
climdex.pcic.test.n.or.x.dates.for.winter.season <- function() {
  test.indices <- names(climdex.min.max.idx.list)[grepl("t", names(climdex.min.max.idx.list))]
  cal <- 365
  test.dates <- seq(as.PCICt("1960-12-01", cal = cal), as.PCICt("1961-02-28", cal = cal), by = "days")
  
  test.tmax <- c(rep(2,2),rep(20, 29), rep(-10, 59))  # Max temperature spikes in December
  test.tmin <- c(rep(2,2),rep(-10, 29), rep(20, 59))  # Min temperature spikes in December
  ci.nx.eoy <- climdexInput.raw(tmax = test.tmax, tmin = test.tmin, tmax.dates = test.dates, tmin.dates = test.dates, base.range=c(1960, 1961))
  
  for (idx in test.indices) {
    data <- get.data.for.idx(ci.nx.eoy, idx)
    
    freq <- "seasonal"
    fun <- paste("climdex", idx, sep = ".")
    
    result <- do.call(fun, list(ci.nx.eoy, freq = freq, include.exact.dates = TRUE))
    expected <- get.n.or.x.result(idx, ci.nx.eoy, freq)
    result$ymd <- as.character(result$ymd)
    checkIdentical(length(expected), nrow(result), paste("Lengths differ. Expected:", length(expected), "Result:", nrow(result)))
    are.not.all.na(expected, result$ymd)
    for (i in seq_along(expected)) {
      expected.val <- data[ci.nx.eoy@dates == expected[[i]]]
      expected.val <- ifelse(length(expected.val) == 0, NA, expected.val)
      checkIdentical(as.character(expected[[i]]), result$ymd[i], paste("Idx:", idx, "Expected:", as.character(expected[[i]]), "Result:", as.character(result$ymd[i])))
      checkTrue(is.almost.equal(as.numeric(expected.val), as.numeric(result$val[i])),
                msg = paste("Idx:", idx, "Expected:", as.numeric(expected.val), "Result:", as.numeric(result$val[i])))
    }

  }
}


# Generic function to get results for rx1day and rx5day indices.
get.Rxnday.result <- function(idx, ci.csv, freq = c("monthly", "annual", "seasonal"), ndays, center.mean.on.last.day) {
  data <- ci.csv@data$prec
  na.mask <- ci.csv@namasks[[freq]]$prec
  fun <- which.max
  if (ndays == 5) {
    data[is.na(data)] <- 0
    prec.runsum <- climdex.pcic:::running.mean(data, ndays)
    prec.runsum[is.na(prec.runsum)] <- 0

    if (center.mean.on.last.day) {
      k2 <- ndays %/% 2
      prec.runsum <- c(rep(0, k2), prec.runsum[1:(length(prec.runsum) - k2)])
    }
    data <- prec.runsum
  }
  factor.extremes <- tapply(data, ci.csv@date.factors[[match.arg(freq)]], fun)
  date.factors <- unique(ci.csv@date.factors[[match.arg(freq)]])
  return(expected.exact.date(ci.csv, data, factor.extremes, date.factors, freq, na.mask))
}


# Test exact dates returned for Rx1day and Rx5day indices.
climdex.pcic.test.exact.date.rxnd.indices <- function() {
  test.indices <- names(climdex.min.max.idx.list)[grepl("r", names(climdex.min.max.idx.list))]
  date.factors <- c("annual", "monthly", "seasonal")

  for (idx in test.indices) {
    ndays <- as.numeric(substr(idx, 3, 3))
    fun <- paste("climdex", idx, sep = ".")

    for (freq in date.factors) {
      result <- do.call(fun, list(ci.csv, freq = freq, include.exact.dates = TRUE))
      center.mean.on.last.day <- FALSE
      expected <- get.Rxnday.result(idx, ci.csv, freq, ndays, center.mean.on.last.day)
      result$ymd <- as.character(result$ymd)
      checkIdentical(length(expected), nrow(result), paste("Lengths Differ. Expected:", length(expected),"Result:", nrow(result)))
      are.not.all.na(expected, result$ymd)
      for (i in seq_along(result$ymd)) {
        if (!is.na(result$ymd[i])) {
          if (ndays == 5) {
            window.start <- expected[[i]] - 2 * 86400
            window.end <- expected[[i]] + 2 * 86400
            expected.val <- sum(ci.csv@data$prec[ci.csv@dates >= window.start & ci.csv@dates <= window.end], na.rm = TRUE)
          } else {
            expected.val <- ci.csv@data$prec[ci.csv@dates == expected[[i]]]
          }
        } else {
          expected.val <- NA
          checkTrue(is.na(expected[[i]]) && is.na(result$val[i]))
        }


        checkIdentical(as.character(expected[[i]]), result$ymd[i], paste("Idx:", idx, "Expected: ", as.character(expected[[i]]), "Result: ", as.character(result$ymd[i])))
        checkTrue(is.almost.equal(as.numeric(expected.val), as.numeric(result$val[i])), 
                  msg = paste("Idx:", idx, "Expected: ", as.numeric(expected.val), "Result: ", as.numeric(result$val[i])))
      }
    }
  }
}
# Check that rx5day works with the mean centered on the last day of the window.
climdex.pcic.test.rx5d.center.mean.on.last.day <- function() {
  date.factors <- c("annual", "monthly", "seasonal")
  ndays <- 5 
  idx <-"rx5day"
  fun <- "climdex.rx5day"
  for (freq in date.factors) {
    center.mean.on.last.day <- TRUE
    result <- do.call(fun, list(ci.csv, freq = freq, center.mean.on.last.day = center.mean.on.last.day, include.exact.dates = TRUE))
    expected <- get.Rxnday.result(idx, ci.csv, freq, ndays, center.mean.on.last.day)
    result$ymd <- as.character(result$ymd)
    checkIdentical(length(expected), nrow(result), paste("Lengths differ. Expected:", length(expected),"Result:", nrow(result)))
    are.not.all.na(expected, result$ymd)
    for (i in seq_along(result$ymd)) {
      if (!is.na(result$ymd[i])) {
        if (ndays == 5) {
          if(center.mean.on.last.day){
            window.start <- expected[[i]] - 4 * 86400
            window.end <- expected[[i]] 
            expected.val <- sum(ci.csv@data$prec[ci.csv@dates >= window.start & ci.csv@dates <= window.end], na.rm = TRUE)
          }
          else{
            window.start <- expected[[i]] - 2 * 86400
            window.end <- expected[[i]] + 2 * 86400
            expected.val <- sum(ci.csv@data$prec[ci.csv@dates >= window.start & ci.csv@dates <= window.end], na.rm = TRUE)
          }

        }
      } else {
        expected.val <- NA
        checkTrue(is.na(expected[[i]]) && is.na(result$val[i]))
      }
      
      
      checkIdentical(as.character(expected[[i]]), result$ymd[i], paste("Idx:", idx, "Expected: ", as.character(expected[[i]]), "Result: ", as.character(result$ymd[i])))
      checkTrue(is.almost.equal(as.numeric(expected.val), as.numeric(result$val[i])), 
                msg = paste("Idx:", idx, "Expected: ", as.numeric(expected.val), "Result: ", as.numeric(result$val[i])))
    }
  }
}

climdex.pcic.test.rx5d.center.mean.on.last.day <- function() {
  date.factors <- c("annual", "monthly", "seasonal")
  ndays <- 5 
  idx <-"rx5day"
  fun <- "climdex.rx5day"
  for (freq in date.factors) {
    center.mean.on.last.day <- TRUE
    result <- do.call(fun, list(ci.csv, freq = freq, center.mean.on.last.day = center.mean.on.last.day, include.exact.dates = TRUE))
    expected <- get.Rxnday.result(idx, ci.csv, freq, ndays, center.mean.on.last.day)
    result$ymd <- as.character(result$ymd)
    checkIdentical(length(expected), nrow(result), paste("Lengths differ. Expected:", length(expected),"Result:", nrow(result)))
    are.not.all.na(expected, result$ymd)
    for (i in seq_along(result$ymd)) {
      if (!is.na(result$ymd[i])) {
        if (ndays == 5) {
          if(center.mean.on.last.day){
            window.start <- expected[[i]] - 4 * 86400
            window.end <- expected[[i]] 
            expected.val <- sum(ci.csv@data$prec[ci.csv@dates >= window.start & ci.csv@dates <= window.end], na.rm = TRUE)
          }
          else{
            window.start <- expected[[i]] - 2 * 86400
            window.end <- expected[[i]] + 2 * 86400
            expected.val <- sum(ci.csv@data$prec[ci.csv@dates >= window.start & ci.csv@dates <= window.end], na.rm = TRUE)
          }

        }
      } else {
        expected.val <- NA
        checkTrue(is.na(expected[[i]]) && is.na(result$val[i]))
      }
      
      
      checkIdentical(as.character(expected[[i]]), result$ymd[i], paste("Idx:", idx, "Expected: ", as.character(expected[[i]]), "Result: ", as.character(result$ymd[i])))
      checkTrue(is.almost.equal(as.numeric(expected.val), as.numeric(result$val[i])), 
                msg = paste("Idx:", idx, "Expected: ", as.numeric(expected.val), "Result: ", as.numeric(result$val[i])))
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
      spell.bounds <- data.frame(start = format(dates[spell$start], "%Y-%m-%d"), end = format(dates[spell$end], "%Y-%m-%d"), duration = ifelse(is.na(ci@namasks$annual$prec[year]), NA, spell$duration))
    } else {
      spell.bounds <- data.frame(start = NA, end = NA, duration = 0)
    }

    return(spell.bounds)
  })

  return(do.call(rbind, spell.boundary))
}

# Generic to compare the expected and climdex-calculated results for the spell tests.
check.spell.results <- function(expected, result, idx) {
  checkIdentical(nrow(expected), nrow(result), paste("Lengths Differ. Expected:", nrow(expected), "Result:", nrow(result)))
  
  for (i in seq_along(result$start)) {
    if (is.na(expected$duration[i])) {
      expected$start[i] <- NA
      expected$end[i] <- NA
    }
    checkIdentical(as.character(expected$start[i]), result$start[i], paste("Start of spells for index:", idx, "do not agree. Expected:", as.character(expected$start[i]), "Result:", result$start[i]))
    checkIdentical(as.character(expected$end[i]), result$end[i], paste("End of spells for index:", idx, "do not agree. Expected:", as.character(expected$end[i]), "Result:", result$end[i]))
    checkTrue(is.almost.equal(as.numeric(expected$duration[i]), as.numeric(result$duration[i])), 
              msg = paste("Idx:", idx, "Expected:",as.character(expected$start[i]), as.numeric(expected$duration[i]), as.character(expected$end[i]), "Result:", result$start[i], as.numeric(result$duration[i]), result$end[i]))
  }
}


# Test cdd and cwd spells with example data set.
climdex.pcic.test.spell.boundaries <- function() {
  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    if (idx == "cdd") {
      result <- climdex.cdd(ci.csv, spells.can.span.years = F, include.exact.dates = TRUE)
    } else {
      result <- climdex.cwd(ci.csv, spells.can.span.years = F, include.exact.dates = TRUE)
    }

    expected <- get.spell.bounds(ci.csv, idx)
    are.not.all.na(expected$start, result$start)
    are.not.all.na(expected$end, result$end)
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
      result <- climdex.cdd(ci.dry, spells.can.span.years = F, include.exact.dates = TRUE)
      ci <- ci.dry
    } else {
      result <- climdex.cwd(ci.wet, spells.can.span.years = F, include.exact.dates = TRUE)
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
      result <- climdex.cdd(ci.wet, spells.can.span.years = F, include.exact.dates = TRUE)
      ci <- ci.wet
    } else {
      result <- climdex.cwd(ci.dry, spells.can.span.years = F, include.exact.dates = TRUE)
      ci <- ci.dry
    }

    expected <- get.spell.bounds(ci, idx)
    check.spell.results(expected, result, idx)
    
    checkTrue(all(is.na(expected$start)))
    checkTrue(all(is.na(result$start)))
    checkTrue(all(is.na(expected$end)))
    checkTrue(all(is.na(result$end)))

  }
}

# Test for spell of 1 day starting from the second last day of the year.
climdex.pcic.test.end.of.year.spell <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1962-01-01", cal = cal), by = "days")

  test.prec.data.dry <- c(rep(0, length(test.dates) - 2), 1, 1)
  ci.dry <- climdexInput.raw(prec = test.prec.data.dry, prec.dates = test.dates)

  test.prec.data.wet <- c(rep(2, length(test.dates) - 2), 0, 0)
  ci.wet <- climdexInput.raw(prec = test.prec.data.wet, prec.dates = test.dates)

  test.indices <- c("cdd", "cwd")
  for (idx in test.indices) {
    if (idx == "cdd") {
      result <- climdex.cdd(ci.wet, spells.can.span.years = F, include.exact.dates = TRUE)
      ci <- ci.wet
    } else {
      result <- climdex.cwd(ci.dry, spells.can.span.years = F, include.exact.dates = TRUE)
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
      result <- climdex.cdd(ci.ran, spells.can.span.years = F, include.exact.dates = TRUE)
    } else {
      result <- climdex.cwd(ci.ran, spells.can.span.years = F, include.exact.dates = TRUE)
    }
    expected <- get.spell.bounds(ci.ran, idx)
    check.spell.results(expected, result, idx)
  }
}
# Check that the start of the (only) spell for 1962 starts in 1961.
<<<<<<< HEAD
climdex.pcic.test.spells.can.span.years <- function() {
  test.indices <- c("cdd", "cwd")
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1962-12-31", cal = cal), by = "days")

  for (idx in test.indices) {
    if (idx == "cdd") {
      test.prec.data <- rep(2, length(test.dates))
      test.prec.data[(cal - 20):(cal - 10)] <- 0
      test.prec.data[(cal - 1):(cal + 4)] <- 0
      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)

      expected <- data.frame(start <- c("1961-12-11", "1961-12-30"), duration <- c(11, 6), end <- c("1961-12-21", "1962-01-04"))
      result <- climdex.cdd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    } else {
      test.prec.data <- rep(0, length(test.dates))
      test.prec.data[(cal - 20):(cal - 10)] <- 2
      test.prec.data[(cal - 1):(cal + 4)] <- 2

      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)

      expected <- data.frame(start <- c("1961-12-11", "1961-12-30"), duration <- c(11, 6), end <- c("1961-12-21", "1962-01-04"))
      result <- climdex.cwd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    }
    check.spell.results(expected, result, idx)
  }
}

# Edge case for different year lengths.
climdex.pcic.test.spells.can.span.leap.year <- function() {
  test.indices <- c("cdd", "cwd")
  cal <- "proleptic_gregorian"
  test.year <- 1964  # Example leap year
  test.dates <- seq(as.PCICt(paste(test.year, "01-01", sep = "-"), cal = cal),
                    as.PCICt(paste(test.year, "12-31", sep = "-"), cal = cal), by = "days")
  cal <- 366
  for (idx in test.indices) {
    if (idx == "cdd") {
      test.prec.data <- rep(2, length(test.dates))
      test.prec.data[40:(cal - 10)] <- 0

      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)
      
      expected <- data.frame(start = c(paste(test.year, "02-09", sep = "-")),
                             duration = c(317),
                             end = c(paste(test.year, "12-21", sep = "-")))
      result <- climdex.cdd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    } else {
      test.prec.data <- rep(0, length(test.dates))
      test.prec.data[40:(cal - 10)] <- 2

      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)
      
      expected <- data.frame(start = c(paste(test.year, "02-09", sep = "-")),
                             duration = c(317),
                             end = c(paste(test.year, "12-21", sep = "-")))
      result <- climdex.cwd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    }
    check.spell.results(expected, result, idx)
  }
}

=======
>>>>>>> 5832ef2 (Add tavg namasks, comments and improve NA seasons test)
climdex.pcic.test.spells.can.span.years <- function() {
  test.indices <- c("cdd", "cwd")
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1962-12-31", cal = cal), by = "days")

  for (idx in test.indices) {
    if (idx == "cdd") {
      test.prec.data <- rep(2, length(test.dates))
      test.prec.data[(cal - 20):(cal - 10)] <- 0
      test.prec.data[(cal - 1):(cal + 4)] <- 0
      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)

      expected <- data.frame(start <- c("1961-12-11", "1961-12-30"), duration <- c(11, 6), end <- c("1961-12-21", "1962-01-04"))
      result <- climdex.cdd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    } else {
      test.prec.data <- rep(0, length(test.dates))
      test.prec.data[(cal - 20):(cal - 10)] <- 2
      test.prec.data[(cal - 1):(cal + 4)] <- 2

      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)

      expected <- data.frame(start <- c("1961-12-11", "1961-12-30"), duration <- c(11, 6), end <- c("1961-12-21", "1962-01-04"))
      result <- climdex.cwd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    }
    check.spell.results(expected, result, idx)
  }
}

# Edge case for different year lengths.
climdex.pcic.test.spells.can.span.leap.year <- function() {
  test.indices <- c("cdd", "cwd")
  cal <- "proleptic_gregorian"
  test.year <- 1964  # Example leap year
  test.dates <- seq(as.PCICt(paste(test.year, "01-01", sep = "-"), cal = cal),
                    as.PCICt(paste(test.year, "12-31", sep = "-"), cal = cal), by = "days")
  cal <- 366
  for (idx in test.indices) {
    if (idx == "cdd") {
      test.prec.data <- rep(2, length(test.dates))
      test.prec.data[40:(cal - 10)] <- 0

      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)
      
      expected <- data.frame(start = c(paste(test.year, "02-09", sep = "-")),
                             duration = c(317),
                             end = c(paste(test.year, "12-21", sep = "-")))
      result <- climdex.cdd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    } else {
      test.prec.data <- rep(0, length(test.dates))
      test.prec.data[40:(cal - 10)] <- 2

      ci <- climdexInput.raw(prec = test.prec.data, prec.dates = test.dates)
      
      expected <- data.frame(start = c(paste(test.year, "02-09", sep = "-")),
                             duration = c(317),
                             end = c(paste(test.year, "12-21", sep = "-")))
      result <- climdex.cwd(ci, spells.can.span.years = TRUE, include.exact.dates = TRUE)
    }
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
expected.gsl <- function(ci, include.exact.dates) {
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
      date.indices <- which(ci@dates >= sy.pcict & ci@dates < ey.pcict)

      tavg <- ci@data$tavg[date.indices]
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
      end.result <- gsl.test.ymd(year, cal, length(tavg)-1, n.h)
      duration <- length(tavg)  - start.idx + 1
    } else {
      end.result <- gsl.test.ymd(year, cal, (ifelse(end.idx > 1, midpoint + end.idx - 1, ifelse(!n.h && next.year.is.leap, midpoint, midpoint + 1))), n.h)
      duration <- difftime(as.Date(end.result),as.Date(start.result))
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

  if (include.exact.dates) {
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
      namask.gsl.monthly <- climdex.pcic:::get.na.mask(gsl.temp.data, gsl.factor.monthly, ci@max.missing.days["annual"])
      ci@namasks$annual$tavg <- climdex.pcic:::get.na.mask(gsl.temp.data, gsl.factor, ci@max.missing.days["annual"]) * as.numeric(tapply(namask.gsl.monthly, gsl.yearmonth.factor, prod))
    }

    ex.r$sl <- ex.r$sl * ci@namasks$annual$tavg
    ex.r$start[is.na(ex.r$sl)] <- NA
    ex.r$end[is.na(ex.r$sl)] <- NA
  }

  return(ex.r)
}
# Generic GSL test function comparing expected with climdex-calculated results.
test.gsl <- function(ci, test.name) {
  expected <- expected.gsl(ci, include.exact.dates = TRUE)
  result <- climdex.gsl(ci, "GSL", include.exact.dates = TRUE)
  checkIdentical(length(expected$start), length(result$start), paste("Lengths differ. Expected:", length(expected$start),"Result:", length(result$start)))
  checkIdentical(length(expected$end), length(result$end), paste("Lengths differ. Expected:", length(expected),"Result:", length(result$end)))
  
  for (i in seq_along(result$start)) {
    if (test.name == "climdex.pcic.test.na.masks.gsl") {
      # Only compare NA rows
      if (i %% 2 == 0) {
        next
      }
    }
    if (test.name != "climdex.pcic.test.no.gsl"){
      are.not.all.na(expected$start, result$start)
      are.not.all.na(expected$end, result$end)
    }
    checkIdentical(expected$start[i], result$start[i], paste("Start of GSL does not match. Expected:", as.character(expected$start[i]), " Result: ", result$start[i], " (", test.name, ")"))
    checkIdentical(expected$end[i], result$end[i], paste("End of GSL does not match. Expected:", as.character(expected$end[i]), " Result: ", result$end[i], " (", test.name, ")"))
    checkTrue(is.almost.equal(as.numeric(expected$sl[i]), as.numeric(result$sl[i])), 
              msg = paste("Idx:GSL\n Year:",rownames(result)[i],"\n Expected: ", as.numeric(expected$sl[i]), "Result: ", as.numeric(result$sl[i])))
  }
}

# Test GSL calculation with example data.
climdex.pcic.test.gsl <- function() {
  test.gsl(ci.csv, "climdex.pcic.test.gsl")
}

# Test for a GSL that spans multiple years.
climdex.pcic.test.multi.year.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
  test.tavg.data <- c(rep(4, cal - 15), rep(6, length(test.dates) - (cal - 15)))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  test.gsl(ci.gsl, "climdex.pcic.test.multi.year.gsl")
}

# Test for GSL that never starts.
climdex.pcic.test.no.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1964-12-31", cal = cal), by = "days")
  test.tavg.data <- c(rep(4, length(test.dates)))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  test.gsl(ci.gsl, "climdex.pcic.test.no.gsl")
}

# Test that the NA masks are applied to results.
climdex.pcic.test.na.masks.gsl <- function() {
  cal <- 365
  test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1969-12-31", cal = cal), by = "days")
  test.dates.factor <- factor(format(test.dates, "%Y"))
  test.tavg.data <- c(sample(1:10, length(test.dates), replace = TRUE))
  ci.gsl <- climdexInput.raw(tavg = test.tavg.data, tavg.dates = test.dates)
  ci.gsl@namasks$annual$tavg <- rep(c(NA, 1), length.out = length(unique(test.dates.factor)))
  test.gsl(ci.gsl, "climdex.pcic.test.na.masks.gsl")
}

# Test GSL in the southern hemisphere.
climdex.pcic.test.gsl.southern.hemisphere <- function() {
  test.gsl(ci.csv.sh, "climdex.pcic.test.gsl.southern.hemisphere")
}

# Get the season from the year and month of an exact date.
format.seasons <- function(months, years) {
  ifelse(months %in% c(12, 1, 2), paste("Winter", as.integer(years) - ifelse(months %in% c(1, 2), 1, 0)),
         ifelse(months %in% 3:5, paste("Spring", years),
                ifelse(months %in% 6:8, paste("Summer", years),
                       ifelse(months %in% 9:11, paste("Fall", years), NA)
                )
         )
  )
}
# Check that each day in results are in the correct year, month or season.
check.single.day.in.factors <- function(freq, result) {
  non.na.rows <- !is.na(result$ymd)
  result.formatted <- switch(
    as.character(freq),
    annual = {format(as.Date(result$ymd[non.na.rows]), "%Y")},
    monthly = {format(as.Date(result$ymd[non.na.rows]), "%Y-%m")},
    seasonal = {
      months <- as.integer(format(as.Date(result$ymd[non.na.rows]), "%m"))
      years <- format(as.Date(result$ymd[non.na.rows]), "%Y")
      format.seasons(months, years)
    }
  )
  checkEquals(result.formatted, rownames(result)[non.na.rows], c("Found the following mismatches:",result.formatted[result.formatted != rownames(result)[non.na.rows]]))
}


# Check that all the bounds of a spell, when spells cannot span years are contained within the bounds of a date factor. 
check.duration.bounds.in.factors <- function(freq, result, spells.can.span.years = F, is.gsl = F) {
  non.na.starts <- !is.na(result$start)
  non.na.ends <- !is.na(result$end)
  
  starts.formatted <- {format(as.Date(result$start[non.na.starts]), "%Y")}
  ends.formatted <- {format(as.Date(result$end[non.na.ends]), "%Y")}
  
  checkEquals(starts.formatted, rownames(result)[non.na.starts], c("Found the following mismatches:",starts.formatted[starts.formatted != rownames(result)[non.na.starts]]))
  checkEquals(ends.formatted, rownames(result)[non.na.ends], c("Found the following mismatches:",ends.formatted[ends.formatted != rownames(result)[non.na.ends]]))

}

# Test that the exact dates for all indices except for GSL have exact dates returned are in their associated date factor.
climdex.pcic.tests.exact.dates.are.in.factors <- function() {
  test.indices <- names(climdex.min.max.idx.list)
  for(idx in test.indices){
    fun <- paste("climdex", idx, sep = ".")
    date.factors <- c("annual", "monthly", "seasonal")
    for (freq in date.factors) {
      result <- do.call(fun, list(ci.csv, freq = freq, include.exact.dates = TRUE))
      check.single.day.in.factors(freq, result)
    }
  }
  
  freq<- "annual"
  result <- climdex.cdd(ci.csv, spells.can.span.years = F, include.exact.dates = TRUE)
  check.duration.bounds.in.factors(freq, result)
  result <- climdex.cwd(ci.csv, spells.can.span.years = F, include.exact.dates = TRUE)
  check.duration.bounds.in.factors(freq, result)

}


checkTypes <- function(result.e.d.vals, result.n.d) {
  for (i in seq_along(result.e.d.vals)){
    checkIdentical(typeof(result.e.d.vals[i]),typeof(result.n.d[i]), paste("Different types: With exact dates:", typeof(result.e.d.vals[i]), "Without:", typeof(result.n.d[i])))
  }
  
}
climdex.pcic.test.consistent.indices.return.types <- function() {
  
  start_date <- as.PCICt("1961-01-01", cal = "365")
  end_date <- as.PCICt("1967-12-30", cal = "365")
  dates <- seq(start_date, end_date, by = "days")
  
  set.seed(123)
  n <- length(dates)
  tmax <- runif(n, -10, 35)
  tmin <- runif(n, -15, 30)
  prec <- runif(n, 0, 40)
  
  na_indices <- sample(1:n, size = round(n * 0.1))
  tmax[na_indices] <- NA
  tmin[na_indices] <- NA
  prec[na_indices] <- NA
  
  # Create climdexInput object
  ci.types.test <- climdexInput.raw(tmax = tmax, tmin = tmin, prec = prec, tmax.dates = dates, tmin.dates = dates, prec.dates = dates)
  test.indices <- c("txx", "tnn", "tnx", "txn", "rx1day", "rx5day")
  for(idx in test.indices){
    fun <- paste("climdex", idx, sep = ".")
    date.factors <- c("annual", "monthly", "seasonal")
    for (freq in date.factors) {
      result.e.d <- do.call(fun, list(ci.types.test, freq = freq, include.exact.dates = TRUE))
      result.n.d <- do.call(fun, list(ci.types.test, freq = freq, include.exact.dates = FALSE))
      checkTypes(result.e.d$val, result.n.d)
    }
  }

  freq<- "annual"
  result.e.d <- climdex.cdd(ci.types.test, spells.can.span.years = F, include.exact.dates = TRUE)
  result.n.d <- climdex.cdd(ci.types.test, spells.can.span.years = F, include.exact.dates = F)
  checkTypes(result.e.d$duration, result.n.d)
  result.e.d <- climdex.cwd(ci.types.test, spells.can.span.years = F, include.exact.dates = TRUE)
  result.n.d <- climdex.cwd(ci.types.test, spells.can.span.years = F, include.exact.dates = F)
  checkTypes(result.e.d$duration, result.n.d)
  result.e.d <- climdex.gsl(ci.types.test,"GSL", include.exact.dates = TRUE)
  result.n.d <- climdex.gsl(ci.types.test,"GSL", include.exact.dates = F)
  checkTypes(result.e.d$sl, result.n.d)
}

# climdex.pcic.test.indicies.on.random.datasets <- function(){
#   for (i in 1:5000){
#     cal <- 365
#     test.dates <- seq(as.PCICt("1961-01-01", cal = cal), as.PCICt("1967-12-31", cal = cal), by = "days")
#     test.prec.ran <- c(sample(0:2, length(test.dates), replace = TRUE))
#     test.tmin.ran <- c(sample(-10:32, length(test.dates), replace = TRUE))
#     test.tmax.ran <- c(sample(-10:32, length(test.dates), replace = TRUE))
#     
# 
#     years <- format(test.dates, "%Y")
#     unique.years <- unique(years)
#     for (year in unique.years) {
#       year.indices <- which(years == year)
#       na.indices <- sample(year.indices, sample(1:7, 1))
#       test.prec.ran[na.indices] <- NA
#       test.tmin.ran[na.indices] <- NA
#       test.tmax.ran[na.indices] <- NA
#     }
#     
#     
#     ci.ran <- climdexInput.raw(tmin =test.tmin.ran, tmax = test.tmax.ran, prec = test.prec.ran, tmin.dates=test.dates, tmax.dates = test.dates, prec.dates = test.dates)
#     
#     test.indices <- c("txx", "tnn", "tnx", "txn", "rx1day", "rx5day")
#     for(idx in test.indices){
#       fun <- paste("climdex", idx, sep = ".")
#       date.factors <- c("annual", "monthly", "seasonal")
#       for (freq in date.factors) {
#         do.call(fun, list(ci.ran, freq = freq, include.exact.dates = TRUE))
#         
#       }
#     }
#     
#     freq<- "annual"
#     climdex.cdd(ci.ran, spells.can.span.years = F, include.exact.dates = TRUE)
#     climdex.cwd(ci.ran, spells.can.span.years = F, include.exact.dates = TRUE)
#     climdex.gsl(ci.ran,"GSL", include.exact.dates = TRUE)
#   } 
#   
#   
# }


checkTypes <- function(result.e.d.vals, result.n.d) {
  for (i in seq_along(result.e.d.vals)){
    checkIdentical(typeof(result.e.d.vals[i]),typeof(result.n.d[i]), paste("Different types: With exact dates:", typeof(result.e.d.vals[i]), "Without:", typeof(result.n.d[i])))
  }
  
}
climdex.pcic.test.consistent.indices.return.types <- function() {
  
  start_date <- as.PCICt("1961-01-01", cal = "365")
  end_date <- as.PCICt("1967-12-30", cal = "365")
  dates <- seq(start_date, end_date, by = "days")
  
  set.seed(123)
  n <- length(dates)
  tmax <- runif(n, -10, 35)
  tmin <- runif(n, -15, 30)
  prec <- runif(n, 0, 40)
  
  na_indices <- sample(1:n, size = round(n * 0.1))
  tmax[na_indices] <- NA
  tmin[na_indices] <- NA
  prec[na_indices] <- NA
  
  # Create climdexInput object
  ci.types.test <- climdexInput.raw(tmax = tmax, tmin = tmin, prec = prec, tmax.dates = dates, tmin.dates = dates, prec.dates = dates)
  test.indices <- names(climdex.min.max.idx.list)
  for(idx in test.indices){
    fun <- paste("climdex", idx, sep = ".")
    date.factors <- c("annual", "monthly", "seasonal")
    for (freq in date.factors) {
      result.e.d <- do.call(fun, list(ci.types.test, freq = freq, include.exact.dates = TRUE))
      result.n.d <- do.call(fun, list(ci.types.test, freq = freq, include.exact.dates = FALSE))
      checkTypes(result.e.d$val, result.n.d)
    }
  }

  freq<- "annual"
  result.e.d <- climdex.cdd(ci.types.test, spells.can.span.years = F, include.exact.dates = TRUE)
  result.n.d <- climdex.cdd(ci.types.test, spells.can.span.years = F, include.exact.dates = F)
  checkTypes(result.e.d$duration, result.n.d)
  result.e.d <- climdex.cwd(ci.types.test, spells.can.span.years = F, include.exact.dates = TRUE)
  result.n.d <- climdex.cwd(ci.types.test, spells.can.span.years = F, include.exact.dates = F)
  checkTypes(result.e.d$duration, result.n.d)
  result.e.d <- climdex.gsl(ci.types.test,"GSL", include.exact.dates = TRUE)
  result.n.d <- climdex.gsl(ci.types.test,"GSL", include.exact.dates = F)
  checkTypes(result.e.d$sl, result.n.d)
}



