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
  var.list <- c("tmax", "tmin", "prec", "tavg")
  for (var in var.list){
    # Test 1.1: Check if the first season is NA based on the start of the time series
    # - If the first month of the date series is not a multiple of 3, expect the first season to be NA
    checkTrue((as.integer(format(ci.csv@dates[1], "%m")) %% 3 != 0) && is.na(ci.csv@namasks$seasonal[[var]][1]))
    # Test 1.2: Check if the last season is NA based on the end of the time series
    # - If the last month of the date series + 1 is not a multiple of 3, expect the last season to be NA
    checkTrue((as.integer(format(ci.csv@dates[length(ci.csv@dates)], "%m")) + 1 %% 3 != 0) && is.na(ci.csv@namasks$seasonal[[var]][length(ci.csv@namasks$season[[var]])]))
    # Test 1.3: Check for NA seasons where conditions in 1.1 and 1.2 do not apply
    # A season is marked as NA if it contains at least one month with missing data or if the sum of missing days within the season is greater than or equal to max.missing.days$season.
    unique_seasons <- unique(ci.csv@date.factors$seasonal)
    unique_months <- unique(ci.csv@date.factors$monthly)
    na_seasons <- unique_seasons[is.na(ci.csv@namasks$seasonal[[var]])]
    na_months <- unique_months[is.na(ci.csv@namasks$monthly[[var]])]
    results <- logical(length(na_seasons))
    
    for (i in seq_along(na_seasons)) {
      season <- na_seasons[i]
      months_of_season <- unique(ci.csv@date.factors$monthly[ci.csv@date.factors$season == season])
      any_month_na <- any(months_of_season %in% na_months)
      results[i] <- any_month_na
      na_days <- sum(is.na(ci.csv@data[[var]][ci.csv@date.factors$season == season]))
      results[i] <- any_month_na || (na_days >= ci.csv@max.missing.days[["seasonal"]])
      
    }
    checkTrue(all(results))
  }
  
}


climdex.pcic.test.seasonal.indices <- function() {
  # Test 2:
  # This section handles various scenarios based on the type of index being calculated:
  # - For non-averaged stats (e.g., min, max), the seasonal value is the min/max of the constituent months.
  # - For averaged stats, the seasonal value may differ from the average of the 3 constituent months due to varying month lengths, so it's compared against a calculated dtr (daily temperature range).
  # - For percentile calculations, each season's percentile is computed against a weighted monthly percentile.

  min_max <- list(
    "tnn" = list(stat = "min", var = "tmin"), "tnx" = list(stat = "max", var = "tmin"),
    "txn" = list(stat = "min", var = "tmax"), "txx" = list(stat = "max", var = "tmax"),
    "rx1day" = list(stat = "max", var = "prec"), "rx5day" = list(stat = "max", var = "prec")
  )
  mean_indices <- list("dtr" = list(var = c("tmin", "tmax")))
  bootstrap_indices <- list("tn10p" = list(var = "tmin"), "tn90p" = list(var = "tmin"), "tx10p" = list(var = "tmax"), "tx90p" = list(var = "tmax"))
  seasonal_indices <- c(names(min_max), names(mean_indices), names(bootstrap_indices))
  ci_subset <- function(ci.csv, dates, season) {
    return(climdexInput.raw(
      tmax = ci.csv@data$tmax[ci.csv@date.factors$seasonal == season],
      tmin = ci.csv@data$tmin[ci.csv@date.factors$seasonal == season],
      prec = ci.csv@data$prec[ci.csv@date.factors$seasonal == season],
      tavg = ci.csv@data$tmin[ci.csv@date.factors$seasonal == season],
      tmax.dates = dates, tmin.dates = dates, prec.dates = dates, tavg.dates = dates, base.range = c(1981, 1990)
    ))
  }
  for (index in seasonal_indices) {
    if (index %in% names(min_max)) {
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
    } else if (index %in% names(mean_indices)) {
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
    } else if (index %in% names(bootstrap_indices)) {
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
    }
  }
}
