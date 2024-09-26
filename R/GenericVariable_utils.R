# Utility function to validate arguments for scalar and vector data.
check.generic.argument.validity <- function( data, dates, max.missing.days) {
  
  stopifnot(length(max.missing.days) == 3 && all(c("annual", "monthly", "seasonal") %in% names(max.missing.days)))
  
  if (!is.numeric(data)) {
    stop("Primary Data must be numeric.")
  }
 
  if (length(data) != length(dates)) {
    stop("Primary data and dates must have the same length.")
  }
  
  if(!is.null(dates) && !inherits(dates, "PCICt"))
      stop(paste("Dates must be of class PCICt."))
}

# Utility function to handle date ranges and generate date factors.
date_info <- function(dates) {
  cal <- attr(dates, "cal")
  
  last.day.of.year <- get.last.monthday.of.year(dates)
  
  date.range <- as.PCICt(paste(as.numeric(format(range(dates), "%Y", tz = "GMT")), c("01-01", last.day.of.year), sep = "-"), cal = cal)
  date.series <- seq(date.range[1], date.range[2], by = "day")
  
  jdays <- get.jdays.replaced.feb29(get.jdays(date.series))
  
  season_with_year <- classify_meteorological_season_with_year(date.series)
  
  date.factors <- list(
    annual = factor(format(date.series, format = "%Y", tz = "GMT")),
    monthly = factor(format(date.series, format = "%Y-%m", tz = "GMT")),
    seasonal = factor(season_with_year, levels = unique(season_with_year))
  )
  
  return(list(
    cal = cal,
    date.series = date.series,
    date.factors = date.factors,
    jdays = jdays
  ))
}

# Generates NA masks based on filled data and date factors
generate_namasks <- function(filled.list, date.factors, max.missing.days) {
  namasks <- list(
    annual = lapply(filled.list, get.na.mask, date.factors$annual, max.missing.days["annual"]),
    monthly = lapply(filled.list, get.na.mask, date.factors$monthly, max.missing.days["monthly"]),
    seasonal = lapply(filled.list, get.na.mask, date.factors$seasonal, max.missing.days["seasonal"]))
    # Vectors: Combine the masks for magnitude and direction
  if ("primary" %in% names(filled.list) && "secondary" %in% names(filled.list)) {
    # Synchronize annual masks
    namasks$annual$primary <- namasks$annual$primary * namasks$annual$secondary
    namasks$annual$secondary <- namasks$annual$primary
    
    # Synchronize monthly masks
    namasks$monthly$primary <- namasks$monthly$primary * namasks$monthly$secondary
    namasks$monthly$secondary <- namasks$monthly$primary
    
    # Synchronize seasonal masks
    namasks$seasonal$primary <- namasks$seasonal$primary * namasks$seasonal$secondary
    namasks$seasonal$secondary <- namasks$seasonal$primary
  }
  namasks$annual <- lapply(names(namasks$annual), function(v) {
    d <- namasks$annual[[v]] * as.numeric(tapply(namasks$monthly[[v]], rep(seq_along(namasks$annual[[v]]), each = 12), prod))
    dimnames(d) <- dim(d) <- NULL
    d
  })
  names(namasks$annual) <- names(namasks$seasonal) <- names(namasks$monthly)
  
  
  season_month_counts <- sapply(unique(date.factors$seasonal), function(season) {
    length(unique(date.factors$monthly[date.factors$seasonal == season]))
  })
  data.vars <- names(filled.list)

  for (var in data.vars) {
    seasonal_namasks <- namasks$seasonal[[var]]
    na_months <- unique(date.factors$monthly)[is.na(namasks$monthly[[var]])]
    seasons_of_na_months <- unique(date.factors$seasonal[date.factors$monthly %in% na_months])
    seasonal_namasks[unique(date.factors$seasonal) %in% seasons_of_na_months] <- NA
    # Identify and set NA for seasons with less than 3 months
    for (season in seq_along(season_month_counts) ) {
      if (!is.na(season_month_counts[season]) && season_month_counts[season] < 3) {
        seasonal_namasks[season] <- NA
      }
    }
    namasks$seasonal[[var]] <- seasonal_namasks
  } 
  return(namasks)
}

generate_filled_list <- function(data, dates, date.series) {
  if (is.vector(data)) {
    return(list(create.filled.series(data, trunc(dates), date.series)))
  } else {
    filled.list <- sapply(data, function(x) { 
      return(create.filled.series(x, trunc(dates), date.series)) 
    }, simplify = FALSE)
    return(filled.list)
  }
}


# Reads data from a CSV file, validates it, and converts date columns to PCICt dates.
read_csv_data <- function(
    file,
    data.columns,
    date.columns,
    date.format,
    na.strings,
    calendar
) {
  # Read the CSV file
  GV.csv <- read.csv(file, na.strings = na.strings)
  
  # Check that data columns exist
  for (col in data.columns) {
    if (!(col %in% names(GV.csv))) {
      stop(paste("Data column", col, "not found in data."))
    }
  }
  
  # Check that date columns exist
  if (!all(date.columns %in% names(GV.csv))) {
    stop(paste("Date columns", paste(date.columns, collapse = ", "), "not found in data."))
  }
  
  # Extract data cols
  data_values <- lapply(data.columns, function(col) GV.csv[[col]])
  
  # Extract the date fields and create date strings
  date_strings <- apply(GV.csv[date.columns], 1, function(row) paste(row, collapse = " "))
  
  # Convert date strings to PCICt dates
  dates <- as.PCICt(strptime(date_strings, format = date.format, tz = "UTC"), cal = calendar)

  return(list(data = data_values, dates = dates))
}
