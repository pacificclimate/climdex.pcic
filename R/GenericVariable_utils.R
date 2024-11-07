# Utility function to validate arguments for scalar and vector data.
check.generic.argument.validity <- function(
    data, 
    dates, 
    max.missing.days, 
    calendar,
    is.vector = FALSE,
    secondary = NULL,
    format = NULL
) {
  # Internal function to validate data and date arguments
  validate_data_dates <- function(data, dates, name) {
    if (missing(data) || is.null(data)) {
      stop(paste(name, "argument is missing."))
    }
    if (length(data) == 0 || length(dates) == 0) {
      stop(paste(name, "and dates must not be empty vectors."))
    }
    if (!is.numeric(data) && (name != "Secondary data")) {
      stop(paste(name, "must be numeric."))
    }
    if (length(data) != length(dates)) {
      stop(paste(name, "and dates must have the same length."))
    }
    if (any(is.na(dates))) {
      stop(paste("Argument 'dates' has NA values."))
    }
  }
  
  # Check max.missing.days
  if (length(max.missing.days) != 3 || !all(c("annual", "monthly", "seasonal") %in% names(max.missing.days))) {
    stop("max.missing.days must be a named vector with 'annual', 'monthly', and 'seasonal' elements.")
  }
  
  # Validate primary data and dates
  validate_data_dates(data, dates, "Primary data")
  
  # Check if dates are PCICt
  if (!inherits(dates, "PCICt")) {
    stop("Dates must be of class PCICt.")
  }
  
  # Vector-specific checks
  if (is.vector) {
    
    validate_data_dates(secondary, dates, "Secondary data")
    # Check that 'format' is provided
    if (missing(format) || is.null(format)) {
      stop("Argument 'format' is missing.")
    }
    
    # Convert the format to lowercase to allow case-insensitive input
    format <- tolower(format)
    
    # Additional validation for format
    if (format %in% c("polar", "cartesian")) {
      if (!is.numeric(secondary)) {
        stop("For 'polar' or 'cartesian' formats, 'secondary' must be numeric.")
      }
    } else if (format == "cardinal") {
      if (!is.character(secondary)) {
        stop("For 'cardinal' format, 'secondary' must be character.")
      }
    } else {
      stop("Invalid 'format'. Use 'polar', 'cartesian', or 'cardinal'.")
    }
  }
  
  # Calendar check: verify it matches one of the recognized types
  valid_calendars <- c("360_day", "360", "365_day", "365", "noleap", "gregorian", "proleptic_gregorian")
  if (!calendar %in% valid_calendars) {
    stop(paste("Invalid calendar type:", calendar, 
               ". Accepted types are '360_day', '360', '365_day', '365', 'noleap', 'gregorian', 'proleptic_gregorian'."))
  }
}

# For single-value-per-month data. Check one day per month and that the day is always the first.
check.single.month.dates <- function(dates) {
  valid_dates <- dates[!is.na(dates)]
  # Check if there is exactly one value per month on the 1st day
  unique_months <- unique(format(valid_dates, "%Y-%m"))
  day_of_month <- as.integer(format(valid_dates, "%d"))
  
  # Check that the length of unique months matches the number of dates, ensuring only one value per month
  if (length(unique_months) != length(valid_dates)) {
    stop("Data must have exactly one value per month.")
  }
  
  # Check that all dates correspond to the 1st day of each month
  if (!all(day_of_month == 1)) {
    stop("Data must be on the 1st day of each month.")
  }
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
  
  calling_func <- as.character(sys.call(-1)[[1]])
  
  # Ensure that the number of data columns matches the type of the calling function
  if (grepl("Scalar", calling_func, ignore.case = TRUE) && length(data.columns) != 1) {
    stop("For scalar data, 'data.columns' should contain exactly 1 column.")
  } else if (grepl("Vector", calling_func, ignore.case = TRUE) && length(data.columns) != 2) {
    stop("For vector data, 'data.columns' should contain exactly 2 columns.")
  }
  
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
