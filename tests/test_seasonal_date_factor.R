library(climdex.pcic)
library(RUnit)
# Test that NA masks work with seasonal aggregation and that indices are correctly calculated.

# Setup
wd <- "./"
ClimVars.tmax <- file.path(wd, "1018935_MAX_TEMP.csv")
ClimVars.tmin <- file.path(wd, "1018935_MIN_TEMP.csv")
ClimVars.prec <- file.path(wd, "1018935_ONE_DAY_PRECIPITATION.csv")

datacol <- list(tmax = "MAX_TEMP", tmin = "MIN_TEMP", prec = "ONE_DAY_PRECIPITATION")
ci.csv <- climdexInput.csv(ClimVars.tmax, ClimVars.tmin, ClimVars.prec, data.columns = datacol, base.range = c(1981, 1990))

unique_seasons <- unique(ci.csv@date.factors$seasonal)
unique_months <- unique(ci.csv@date.factors$monthly)

climdex.pcic.test.seasonal.na.cases <- function() {
  
  start_date <- as.PCICt("1961-01-01", cal = "365")
  end_date <- as.PCICt("1967-12-31", cal = "365")
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
  ci.na.test <- climdexInput.raw(tmax = tmax, tmin = tmin, prec = prec, tmax.dates = dates, tmin.dates = dates, prec.dates = dates)
  ci.na.test
  
  var.list <- c("tmax", "tmin", "prec", "tavg")
  for (var in var.list) {
    # Check if the first season is NA based on the start of the time series
    is_not_start_of_season <- ifelse(!is.na(ci.na.test@dates[1]),
                                     (as.integer(format(ci.na.test@dates[1], "%m")) %% 3 != 0),
                                     TRUE)
    first_season_na <- is.na(ci.na.test@namasks$seasonal[[var]][1])
    checkTrue(is_not_start_of_season == first_season_na, 
              msg = paste("Start of season NA check failed:", is_not_start_of_season, first_season_na))
    
    # Check if the last season is NA based on the end of the time series
    is_not_end_of_season <- ifelse(!is.na(ci.na.test@dates[length(ci.na.test@dates)]),
                                   (as.integer(format(ci.na.test@dates[length(ci.na.test@dates)], "%m"))+1 %% 3 != 0),
                                   TRUE)
    last_season_na <- is.na(ci.na.test@namasks$seasonal[[var]][length(ci.na.test@namasks$season[[var]])])
    checkTrue(is_not_end_of_season == last_season_na, 
              msg = paste("End of season NA check failed:", is_not_end_of_season, last_season_na))
    
    
    # Check for NA seasons where conditions in 1.1 and 1.2 do not apply
    # A season is marked as NA if it contains at least one month with missing data or if the sum of missing days within the season is greater than or equal to max.missing.days$season.
    unique_seasons <- unique(ci.na.test@date.factors$seasonal)
    unique_months <- unique(ci.na.test@date.factors$monthly)
    na_seasons <- unique_seasons[is.na(ci.na.test@namasks$seasonal[[var]])]
    na_months <- unique_months[is.na(ci.na.test@namasks$monthly[[var]])]
    results <- logical(length(na_seasons))

    for (i in seq_along(na_seasons)) {
      season <- na_seasons[i]
      months_of_season <- unique(ci.na.test@date.factors$monthly[ci.na.test@date.factors$season == season])
      any_month_na <- any(months_of_season %in% na_months)
      results[i] <- any_month_na
      na_days <- sum(is.na(ci.na.test@data[[var]][ci.na.test@date.factors$season == season]))
      results[i] <- any_month_na || (na_days >= ci.na.test@max.missing.days[["seasonal"]])
    }
    # Exclude first and last season, as they have been compared in the previous checks
    checkTrue(all(results[-c(1, length(results))]))
  }
}

# Test Seasonal Indices:
# This section handles various scenarios based on the type of index being calculated:
# - For non-averaged stats (e.g., min, max), the seasonal value is the min/max of the constituent months.
# - For averaged stats, the seasonal value may differ from the average of the 3 constituent months due to varying month lengths, so it's compared against a calculated dtr (daily temperature range).
# - For percentile calculations, each season's percentile is computed against a weighted monthly percentile.

ci_subset <- function(ci.csv, dates, season) {
  return(climdexInput.raw(
    tmax = ci.csv@data$tmax[ci.csv@date.factors$seasonal == season],
    tmin = ci.csv@data$tmin[ci.csv@date.factors$seasonal == season],
    prec = ci.csv@data$prec[ci.csv@date.factors$seasonal == season],
    tavg = ci.csv@data$tmin[ci.csv@date.factors$seasonal == season],
    tmax.dates = dates, tmin.dates = dates, prec.dates = dates, tavg.dates = dates, base.range = c(1981, 1990)
  ))
}

climdex.pcic.test.seasonal.min.max.indices <- function() {
  min_max <- climdex.min.max.idx.list
  lapply(names(min_max), function(index) {
    stat <- min_max[[index]]$stat
    varcol <- min_max[[index]]$var
    data_seasons <- unique_seasons[!is.na(ci.csv@namasks$seasonal[[varcol]])]

    for (season in data_seasons) {
      dates <- unique(ci.csv@dates[ci.csv@date.factors$seasonal == season])
      ci_season <- ci_subset(ci.csv, dates, season)
      fun <- paste("climdex", index, sep = ".")
      clim.result <- get(fun)(ci_season, freq = "seasonal")
      if (index == "rx5day") {
        precip <- (ci_season@data[[varcol]][ci_season@date.factors$seasonal == season])
        statval <- get(stat)(sapply(1:(length(precip) - 4), function(i) sum(precip[i:(i + 4)])), na.rm = TRUE)
      } else {
        statval <- get(stat)(ci_season@data[[varcol]][ci_season@date.factors$seasonal == season], na.rm = TRUE)
      }
      checkEqualsNumeric(statval, clim.result[season])
    }
  })
}

climdex.pcic.test.seasonal.mean.indices <- function() {
  mean_indices <- climdex.mean.idx.list
  index <- "dtr"
  vars <- mean_indices[[index]]$var
  dtr_vars <- pmin(!is.na(ci.csv@namasks$seasonal[[vars[1]]]), !is.na(ci.csv@namasks$seasonal[[vars[2]]]), na.rm = TRUE)
  data_seasons <- unique_seasons[dtr_vars > 0]
  for (season in data_seasons) {
    dates <- unique(ci.csv@dates[ci.csv@date.factors$seasonal == season])
    ci_season <- ci_subset(ci.csv, dates, season)
    fun <- paste("climdex", index, sep = ".")
    clim.result <- get(fun)(ci_season, freq = "seasonal")
    dtrl <- na.omit((ci.csv@data$tmax[ci.csv@date.factors$seasonal == season] - ci.csv@data$tmin[ci.csv@date.factors$seasonal == season]))
    dtrm <- sum(as.numeric(dtrl)) / length(as.numeric(dtrl))
    checkEqualsNumeric(clim.result[season], dtrm, tolerance = 1e-6)
  }
}

climdex.pcic.test.seasonal.bootstrap.indices <- function() {
  bootstrap_indices <- climdex.bootstrap.idx.list
  lapply(names(bootstrap_indices), function(index) {
    fun <- paste("climdex", index, sep = ".")
    varcol <- bootstrap_indices[[index]]$var
    data_seasons <- unique_seasons[!is.na(ci.csv@namasks$seasonal[[varcol]])]
    clim.result <- get(fun)(ci.csv, freq = "seasonal")
    clim.result.monthly <- get(fun)(ci.csv, freq = "monthly")
    for (season in data_seasons) {
      months_of_season <- unique(ci.csv@date.factors$monthly[ci.csv@date.factors$seasonal == season])
      days_in_months <- sapply(months_of_season, function(month) sum(!is.na(ci.csv@data[[varcol]][ci.csv@date.factors$monthly == month])))
      monthly_weight <- sapply(days_in_months, function(month) month / sum(days_in_months))
      monthly_percentile <- sum(clim.result.monthly[months_of_season] * monthly_weight)
      checkEqualsNumeric(clim.result[season], monthly_percentile, tolerance = 1e-1)
    }
  })
}