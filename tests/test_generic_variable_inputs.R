climdex.pcic.test.na_masks_with_missing_data_thresholds <- function() climdex.pcic.test.na_masks_with_missing_data_thresholds <- function() {
  set.seed(123)
  speed <- runif(366, 0, 20)
  direction <- runif(366, 0, 360)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 366)
  
  # Introduce missing data exceeding monthly and seasonal thresholds
  speed[1:10] <- NA        # First 10 days of January missing
  speed[61:81] <- NA       # 21 days missing in March
  speed[183:203] <- NA     # 21 days missing in July
  
  # Define max missing days thresholds
  max_missing_days <- c(annual = 15, monthly = 5, seasonal = 15)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = max_missing_days,
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  # Function to calculate number of missing days per period
  calculate_missing_days <- function(data, date.factor) {
    tapply(is.na(data), date.factor, sum)
  }
  
  # Calculate missing days for primary data
  monthly_missing_days_primary <- calculate_missing_days(vector_obj@primary, vector_obj@date.factors$monthly)
  seasonal_missing_days_primary <- calculate_missing_days(vector_obj@primary, vector_obj@date.factors$seasonal)
  annual_missing_days_primary <- sum(is.na(vector_obj@primary))
  
  # Expected months and seasons to be NA due to missing data
  months_expected_na <- names(monthly_missing_days_primary[monthly_missing_days_primary > max_missing_days["monthly"]])
  seasons_expected_na <- names(seasonal_missing_days_primary[seasonal_missing_days_primary > max_missing_days["seasonal"]])
  
  # Adjust for the extra season (Winter-2019) that is automatically marked as NA
  all_seasons <- levels(vector_obj@date.factors$seasonal)
  first_season <- all_seasons[1]
  last_season <- all_seasons[length(all_seasons)]
  
  # Combine expected NA seasons from missing data and incomplete seasons
  total_expected_na_seasons <- unique(c(first_season,seasons_expected_na, last_season))
  
  # Check that the total number of NA seasons matches the expected count
  seasonal_namasks_primary <- vector_obj@namasks$seasonal$primary
  actual_na_seasons <- all_seasons[is.na(seasonal_namasks_primary)]
  
  checkEquals(length(actual_na_seasons), length(total_expected_na_seasons),
              msg = paste("Total number of NA seasons does not match expected (expected ", length(total_expected_na_seasons), ").", sep=""))
  
  # Check that actual and expected NA seasons are identical
  checkTrue(identical(actual_na_seasons,total_expected_na_seasons), msg = "Actual NA seasons do not match expected")
  
  # Check that the specific seasons are marked NA
  checkTrue(all(is.na(seasonal_namasks_primary[total_expected_na_seasons])),
            msg = "Seasons expected to be NA are not marked NA.")
  
  # Check that the annual namask is NA
  annual_namask_primary <- vector_obj@namasks$annual$primary
  checkTrue(is.na(annual_namask_primary),
            msg = "Annual namask is not NA when missing data exceeds annual threshold.")
  
  # Verify total counts of NA months
  monthly_namasks_primary <- vector_obj@namasks$monthly$primary
  
  checkEquals(sum(is.na(monthly_namasks_primary)), length(months_expected_na),
              msg = paste("Total number of NA months does not match expected (", length(months_expected_na), " months).", sep=""))
  
  # Check that the specific months are marked NA
  checkTrue(all(is.na(monthly_namasks_primary[months_expected_na])),
            msg = "Months expected to be NA are not marked NA.")
}



climdex.pcic.test.scalar.raw.and.csv.construction <- function() {
  set.seed(123)

  # Raw
  data <- runif(366, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 366)

  scalar_obj_raw <- climdexGenericScalar.raw(
    data = data,
    dates =dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  checkEquals(length(scalar_obj_raw@data), length(scalar_obj_raw@dates), "Raw scalar construction data length does not match dates.")

  # CSV
  csv_data <- data.frame(date = as.character(dates), data = data)
  temp_csv <- tempfile()
  write.csv(csv_data, temp_csv, row.names = FALSE)

  scalar_obj_csv <- climdexGenericScalar.csv(
    file = temp_csv,
    data.column = "data",
    date.columns = "date",
    date.format = "%Y-%m-%d",
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  checkEquals(length(scalar_obj_csv@data), length(scalar_obj_csv@dates), "CSV scalar construction data length does not match dates.")
  checkEquals(scalar_obj_raw@dates, scalar_obj_csv@dates, "Date mismatch between raw and CSV scalar objects.")
  checkTrue(all.equal(scalar_obj_csv, scalar_obj_raw), msg = "Scalar_obj built from csv is not identical to raw")
}

climdex.pcic.test.calendar.handling <- function() {
  set.seed(123)
  
  # Create data for two different calendar types
  data_gregorian <- runif(366, 0, 20)
  dates_gregorian <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 366)
  
  data_noleap <- runif(365, 0, 20)
  dates_noleap <- seq(as.PCICt("2020-01-01", cal = "noleap"), by = "day", length.out = 365)
  
  # Gregorian calendar
  scalar_obj_gregorian <- climdexGenericScalar.raw(
    data = data_gregorian,
    dates = dates_gregorian,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # No-leap calendar
  scalar_obj_noleap <- climdexGenericScalar.raw(
    data = data_noleap,
    dates = dates_noleap,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "noleap" # uses a 365 calendar
  )
  
  checkEquals(attr(scalar_obj_gregorian@dates, "cal"), "proleptic_gregorian", msg = "Gregorian calendar mismatch.")
  checkEquals(attr(scalar_obj_noleap@dates, "cal"), "365", msg = "No-leap calendar mismatch.")
  # Verify that February 29th is included in the Gregorian dates
  checkTrue(any(format(scalar_obj_gregorian@dates, "%m-%d") == "02-29"), "Leap day missing in Gregorian calendar dates.")
  # Verify that February 29th is not included in the no-leap dates
  checkTrue(!any(format(scalar_obj_noleap@dates, "%m-%d") == "02-29"), "Leap day present in no-leap calendar dates.")
}

climdex.pcic.test.vector.raw.and.csv.construction <- function() {
  set.seed(123)
  
  # Raw vector data construction
  primary <- runif(366, 0, 20) # Primary component (e.g., magnitude)
  secondary <- runif(366, 0, 360) # Secondary component (e.g., direction)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 366)
  
  vector_obj_raw <- climdexGenericVector.raw(
    primary = primary,
    secondary = secondary,
    dates = dates,
    format = "polar",
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  checkEquals(length(vector_obj_raw@primary), length(vector_obj_raw@dates), "Raw vector construction primary data length does not match dates.")
  checkEquals(length(vector_obj_raw@secondary), length(vector_obj_raw@dates), "Raw vector construction secondary data length does not match dates.")
  
  # CSV vector data construction
  csv_data <- data.frame(date = as.character(dates), primary = primary, secondary = secondary)
  temp_csv <- tempfile()
  write.csv(csv_data, temp_csv, row.names = FALSE)
  
  vector_obj_csv <- climdexGenericVector.csv(
    file = temp_csv,
    primary.column = "primary",
    secondary.column = "secondary",
    date.columns = "date",
    date.format = "%Y-%m-%d",
    format = "polar",
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  checkEquals(length(vector_obj_csv@primary), length(vector_obj_csv@dates), "CSV vector construction primary data length does not match dates.")
  checkEquals(length(vector_obj_csv@secondary), length(vector_obj_csv@dates), "CSV vector construction secondary data length does not match dates.")
  checkEquals(vector_obj_raw@dates, vector_obj_csv@dates, "Date mismatch between raw and CSV vector objects.")
  checkTrue(all.equal(vector_obj_csv, vector_obj_raw), msg = "Vector_obj built from csv is not identical to raw")
}


climdex.pcic.test.vector.construction.validity <- function() {
  set.seed(123)
  
  # Valid data for 'polar' format
  speed <- c(5, 10, 15)
  direction <- c(0, 90, 180)
  dates <- as.PCICt(c("2020-01-01", "2020-01-05", "2020-01-10"), cal = "gregorian")
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Extract the indices of the dates in the filled date series
  date_series <- vector_obj@dates
  indices <- match(dates, date_series)
  
  # Check that data is correctly assigned at the correct positions
  checkEquals(vector_obj@primary[indices], speed, "Primary data not correctly assigned.")
  checkEquals(vector_obj@secondary[indices], direction, "Secondary data not correctly assigned.")
  checkEquals(vector_obj@dates[indices], dates, "Dates not correctly assigned.")
  
  # Check that other positions are NA
  other_indices <- setdiff(seq_along(date_series), indices)
  checkTrue(all(is.na(vector_obj@primary[other_indices])), "Other primary data should be NA.")
  checkTrue(all(is.na(vector_obj@secondary[other_indices])), "Other secondary data should be NA.")
  
  # Check that format is correctly set
  checkEquals(vector_obj@format, "polar", "Format not correctly set.")
  
  # For 'cartesian' format, verify conversion to polar coordinates
  u <- c(5, 0, -5)
  v <- c(0, 5, 0)
  cartesian_dates <- as.PCICt(c("2020-01-02", "2020-01-06", "2020-01-11"), cal = "gregorian")
  
  vector_obj_cartesian <- climdexGenericVector.raw(
    primary = u,
    secondary = v,
    dates = cartesian_dates,
    format = "cartesian",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Convert cartesian to polar manually
  expected_polar <- convert_cartesian_to_polar(u, v)
  
  # Extract indices corresponding to the cartesian dates
  indices_cartesian <- match(cartesian_dates, vector_obj_cartesian@dates)
  
  # Extract the primary and secondary data at those indices
  primary_cartesian <- vector_obj_cartesian@primary[indices_cartesian]
  secondary_cartesian <- vector_obj_cartesian@secondary[indices_cartesian]
  
  # Convert the object's data to polar coordinates
  polar_data <- convert_cartesian_to_polar(primary_cartesian, secondary_cartesian)
  
  # Check that internal conversion matches expected values
  checkEqualsNumeric(polar_data$speed, expected_polar$speed, tolerance = 1e-6, "Magnitude conversion from cartesian to polar incorrect.")
  checkEqualsNumeric(polar_data$direction, expected_polar$direction, tolerance = 1e-6, "Direction conversion from cartesian to polar incorrect.")
}

climdex.pcic.test.scalar_construction_with_malformed_inputs <- function() {
  # Valid data
  data <- c(10.5, 12.3, 11.2)
  dates <- as.PCICt(c("2020-01-01", "2020-01-02", "2020-01-03"), cal = "gregorian")
  
  # Test 1: Non-numeric data
  invalid_data <- c("high", "medium", "low")
  checkException(
    climdexGenericScalar.raw(
      data = invalid_data,
      dates = dates
    ),
    msg = "Error not raised for non-numeric 'data'."
  )
  
  # Test 2: Missing required arguments
  # Missing 'data' argument
  checkException(
    climdexGenericScalar.raw(
      dates = dates
    ),
    msg = "Error not raised when 'data' argument is missing."
  )
  
  # Missing 'dates' argument
  checkException(
    climdexGenericScalar.raw(
      data = data
    ),
    msg = "Error not raised when 'dates' argument is missing."
  )
  
  # Test 3: Mismatched lengths between data and dates
  data_short <- data[1:2]
  checkException(
    climdexGenericScalar.raw(
      data = data_short,
      dates = dates
    ),
    msg = "Error not raised for mismatched lengths between 'data' and 'dates'."
  )
  
  # Test 4: Invalid dates (non-PCICt objects)
  invalid_dates <- as.Date(c("2020-01-01", "2020-01-02", "2020-01-03"))
  checkException(
    climdexGenericScalar.raw(
      data = data,
      dates = invalid_dates
    ),
    msg = "Error not raised for invalid 'dates' data type."
  )
}



climdex.pcic.test.vector_construction_with_malformed_inputs <- function() {
  # Setup valid data for comparison
  valid_speed <- c(5, 10, 15)
  valid_direction_numeric <- c(0, 90, 180)
  valid_direction_cardinal <- c("N", "E", "S")
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 3)
  
  # Test 1: Mismatched lengths between primary and secondary
  speed_short <- c(5, 10)
  checkException(
    climdexGenericVector.raw(
      primary = speed_short,
      secondary = valid_direction_numeric,
      dates = dates,
      format = "polar"
    ),
    msg = "Error not raised for mismatched lengths between primary and secondary."
  )
  
  # Test 2: Invalid format value
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      secondary = valid_direction_numeric,
      dates = dates,
      format = "invalid_format"
    ),
    msg = "Error not raised for invalid format value."
  )
  
  # Test 3: Non-numeric secondary data for 'polar' format
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      secondary = valid_direction_cardinal,
      dates = dates,
      format = "polar"
    ),
    msg = "Error not raised for non-numeric secondary data in 'polar' format."
  )
  
  # Test 4: Non-character secondary data for 'cardinal' format
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      secondary = valid_direction_numeric,
      dates = dates,
      format = "cardinal"
    ),
    msg = "Error not raised for non-character secondary data in 'cardinal' format."
  )
  
  # Test 5: Non-numeric primary data
  invalid_speed <- c("fast", "medium", "slow")
  checkException(
    climdexGenericVector.raw(
      primary = invalid_speed,
      secondary = valid_direction_numeric,
      dates = dates,
      format = "polar"
    ),
    msg = "Error not raised for non-numeric primary data."
  )
  
  # Test 6: Missing required arguments
  # Missing 'secondary' argument
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      dates = dates,
      format = "polar"
    ),
    msg = "Error not raised when 'secondary' argument is missing."
  )
  
  # Missing 'primary' argument
  checkException(
    climdexGenericVector.raw(
      secondary = valid_direction_numeric,
      dates = dates,
      format = "polar"
    ),
    msg = "Error not raised when 'primary' argument is missing."
  )
  
  # Test 7: Mismatched lengths between data and dates
  dates_short <- dates[1:2]
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      secondary = valid_direction_numeric,
      dates = dates_short,
      format = "polar"
    ),
    msg = "Error not raised for mismatched lengths between data and dates."
  )
  
  # Test 8: Invalid dates (non-PCICt objects)
  invalid_dates <- as.Date(c("2020-01-01", "2020-01-02", "2020-01-03"))
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      secondary = valid_direction_numeric,
      dates = invalid_dates,
      format = "polar"
    ),
    msg = "Error not raised for invalid 'dates' data type."
  )
  
  # Test 9: Invalid calendar type
  checkException(
    climdexGenericVector.raw(
      primary = valid_speed,
      secondary = valid_direction_cardinal,
      dates = dates,
      format = "cardinal",
      calendar = "greg"
    ),
    msg = "Error not raised for invalid calendar type."
  )
  
}

