library(climdex.pcic)
library(RUnit)

climdex.pcic.test.single.monthly.scalar.raw.and.csv.construction <- function() {
  set.seed(123)

  scalar_data <- runif(12, 0, 20) # One value per month
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)
  scalar_obj_raw <- climdexSingleMonthlyScalar.raw(
    data = scalar_data,
    dates = dates,
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  csv_data <- data.frame(date = as.character(dates), data = scalar_data)
  temp_csv <- tempfile()
  write.csv(csv_data, temp_csv, row.names = FALSE)

  scalar_obj_csv <- climdexSingleMonthlyScalar.csv(
    file = temp_csv,
    data.column = "data",
    date.columns = "date",
    date.format = "%Y-%m-%d",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  checkEquals(scalar_obj_raw@dates, scalar_obj_csv@dates, "Date mismatch between raw and CSV scalar objects.")
  checkTrue(all.equal(scalar_obj_csv, scalar_obj_raw), msg = "Scalar_obj built from CSV is not identical to raw")
}

climdex.pcic.test.single.monthly.vector.raw.and.csv.construction <- function() {
  set.seed(123)

  primary_data <- runif(12, 0, 20) # One value per month
  secondary_data <- runif(12, 0, 360)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  vector_obj_raw <- climdexSingleMonthlyVector.raw(
    primary = primary_data,
    secondary = secondary_data,
    dates = dates,
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  csv_data <- data.frame(date = as.character(dates), primary = primary_data, secondary = secondary_data)
  temp_csv <- tempfile()
  write.csv(csv_data, temp_csv, row.names = FALSE)

  vector_obj_csv <- climdexSingleMonthlyVector.csv(
    file = temp_csv,
    primary.column = "primary",
    secondary.column = "secondary",
    date.columns = "date",
    date.format = "%Y-%m-%d",
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  checkEquals(vector_obj_raw@dates, vector_obj_csv@dates, "Date mismatch between raw and CSV vector objects.")
  checkTrue(all.equal(vector_obj_csv, vector_obj_raw), msg = "Vector_obj built from CSV is not identical to raw")
}

climdex.pcic.test.SingleMonthlyScalar.raw.missing <- function() {
  set.seed(123)

  # Single monthly value data with an NA value
  data <- c(1:5, NA, 7:12)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  # Create the climdexSingleMonthlyScalar object and expect it to pass without failure
  result <- try(
    scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
}

climdex.pcic.test.SingleMonthlyVector.raw.missing <- function() {
  set.seed(123)

  # Single monthly value vector data with an NA value in the primary component
  primary <- c(1:5, NA, 7:12)
  secondary <- runif(12, 0, 360)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  # Create the climdexSingleMonthlyVector object and expect it to pass without failure
  result <- try(
    vector_obj <- climdexSingleMonthlyVector.raw(primary, secondary, dates, format = "polar", northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.zero.length <- function() {
  # Empty data and dates vectors
  data <- numeric(0)
  dates <- as.PCICt(character(0), cal = "gregorian")

  error_message <- tryCatch(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    error = function(e) e$message
  )

  # Check error message
  checkTrue(
    grepl("Primary data and dates must not be empty vectors.", error_message),
    "Error message is not informative for empty input vectors."
  )
}

climdex.pcic.test.MultiyearScalarContinuous <- function() {
  set.seed(123)

  data <- runif(36, 0, 20) # One value per month for 3 years
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 36)

  # Create the scalar object
  scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian")

  # Validate that dates and data are maintained correctly
  sorted_indices <- order(dates)
  sorted_dates <- dates[sorted_indices]
  sorted_data <- data[sorted_indices]

  obj_dates <- scalar_obj@dates[!is.na(scalar_obj@data)]
  obj_data <- scalar_obj@data[!is.na(scalar_obj@data)]

  checkEquals(obj_data, sorted_data, "Multiyear scalar data is not aligned correctly with sorted dates.")
}

climdex.pcic.test.MultiyearVectorContinuous <- function() {
  set.seed(123)

  primary <- runif(36, 0, 20) # One value per month for 3 years
  secondary <- runif(36, 0, 360)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 36)

  vector_obj <- climdexSingleMonthlyVector.raw(primary, secondary, dates, format = "polar", northern.hemisphere = TRUE, calendar = "gregorian")

  checkEquals(length(vector_obj@primary), length(vector_obj@dates), "Multiyear vector primary length does not match dates.")

  sorted_indices <- order(dates)
  sorted_dates <- dates[sorted_indices]
  sorted_primary <- primary[sorted_indices]
  sorted_secondary <- secondary[sorted_indices]

  obj_dates <- vector_obj@dates[!is.na(vector_obj@primary)]
  obj_primary <- vector_obj@primary[!is.na(vector_obj@primary)]
  obj_secondary <- vector_obj@secondary[!is.na(vector_obj@primary)]

  checkEquals(obj_primary, sorted_primary, "Multiyear vector primary data is not aligned correctly with sorted dates.")
  checkEquals(obj_secondary, sorted_secondary, "Multiyear vector secondary data is not aligned correctly with sorted dates.")
}

climdex.pcic.test.MultiyearWithGaps <- function() {
  set.seed(123)

  # Scalar data for three years with some missing months
  data <- c(runif(11, 0, 20), NA, runif(11, 0, 20), NA, runif(10, 0, 20), NA, NA) # 36 data points with some NA values to align with dates
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 36)
  # Test that the scalar object is built without errors
  result <- try(
    scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
  # Validate non-NA data alignment
  obj_data <- scalar_obj@data[!is.na(scalar_obj@data)]
  obj_dates <- scalar_obj@dates[!is.na(scalar_obj@data)]
  valid_indices <- !is.na(data)
  checkEquals(obj_data, data[valid_indices], "Multiyear scalar data with gaps is not aligned correctly.")
}

climdex.pcic.test.SingleMonthlyScalar.raw.dates.not.first.day <- function() {
  set.seed(123)

  data <- runif(12, 0, 20)
  dates <- seq(as.PCICt("2020-01-02", cal = "gregorian"), by = "month", length.out = 12)

  # Expect an error due to dates not being on the first day
  checkException(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    "Function did not raise an error when dates were not on the first day of the month."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.mismatched.lengths <- function() {
  set.seed(123)

  # Data and dates of different lengths
  data <- runif(11, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  # Expect an error due to mismatched lengths
  checkException(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    "Function did not raise an error when data and dates lengths were mismatched."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.non.numeric.data <- function() {
  data <- rep("non-numeric", 12)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  # Expect an error due to non-numeric data
  checkException(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    "Function did not raise an error when non-numeric data was provided."
  )
}

climdex.pcic.test.SingleMonthlyScalar.csv.invalid.date.format <- function() {
  # Invalid date formats
  data <- runif(12, 0, 20)
  dates <- format(seq(as.Date("2020-01-01"), by = "month", length.out = 12), "%Y/%m/%d") # Different date format
  csv_data <- data.frame(date = dates, data = data)
  temp_csv <- tempfile()
  write.csv(csv_data, temp_csv, row.names = FALSE)

  # Expect an error due to invalid date format
  checkException(
    climdexSingleMonthlyScalar.csv(
      file = temp_csv,
      data.column = "data",
      date.columns = "date",
      date.format = "%Y-%m-%d", # Expecting different format
      northern.hemisphere = TRUE,
      calendar = "gregorian"
    ),
    "Function did not raise an error when invalid date format was provided in CSV."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.invalid.calendar <- function() {
  set.seed(123)

  # Create valid data and dates
  data <- runif(12, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  # Expect an error due to invalid calendar type
  checkException(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "invalid_calendar"),
    "Function did not raise an error when an invalid calendar type was provided."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.NA.dates <- function() {
  set.seed(123)

  # NA in dates
  data <- runif(12, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)
  dates_char <- as.character(dates)
  dates_char[6] <- NA # Insert NA
  dates <- as.PCICt(dates_char, cal = "gregorian")
  # Capture the error message
  error_message <- tryCatch(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    error = function(e) e$message
  )
  error_message
  # Check error message
  checkTrue(
    grepl("Argument 'dates' has NA values.", error_message),
    "Error message is not informative for NA values in dates."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.leap.year <- function() {
  set.seed(123)

  # Leap year
  data <- runif(24, 0, 20)
  dates <- seq(as.PCICt("2019-01-01", cal = "gregorian"), by = "month", length.out = 24)

  result <- try(
    scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
  # Ensure that February 29, 2020, is included
  checkTrue(any(format(scalar_obj@dates, "%Y-%m-%d") == "2020-02-29"), "Leap day not included in dates.")
}

climdex.pcic.test.SingleMonthlyScalar.raw.extreme.values <- function() {
  # Extreme values
  data <- c(-1e10, runif(10, -1e5, 1e5), 1e10)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian")
  scalar_max <- unname(compute.stat.scalar(scalar_obj, "max", "annual", FALSE))
  scalar_min <- unname(compute.stat.scalar(scalar_obj, "min", "annual", FALSE))
  scalar_sum <- unname(compute.stat.scalar(scalar_obj, "sum", "annual", FALSE))
  checkEquals(scalar_max, 1e10, "annual max stat for single monthly scalar was not equal to max of input data")
  checkEquals(scalar_min, -1e10, "annual min stat for single monthly scalar was not equal to min of input data")
  checkEquals(scalar_sum, sum(data), "annual sum stat for single monthly scalar was not equal to sum of input data")
}

climdex.pcic.test.SingleMonthlyVector.raw.missing.secondary <- function() {
  set.seed(123)

  # NA in secondary component
  primary <- runif(12, 0, 20)
  secondary <- c(runif(5, 0, 360), NA, runif(6, 0, 360))
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  result <- try(
    vector_obj <- climdexSingleMonthlyVector.raw(primary, secondary, dates, format = "polar", northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
  checkEquals(vector_obj@primary[!is.na(vector_obj@primary)], vector_obj@primary[!is.na(vector_obj@secondary)])
}


climdex.pcic.test.SingleMonthlyVector.raw.invalid.format <- function() {
  set.seed(123)

  # Valid data
  primary <- runif(12, 0, 20)
  secondary <- runif(12, 0, 360)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  # Expect an error due to invalid format
  checkException(
    climdexSingleMonthlyVector.raw(primary, secondary, dates, format = "invalid_format", northern.hemisphere = TRUE, calendar = "gregorian"),
    "Function did not raise an error when an invalid format was provided."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.error.messages <- function() {
  # Multiple data values per month
  data <- runif(24, 0, 20)
  dates <- c(
    seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12),
    seq(as.PCICt("2020-01-15", cal = "gregorian"), by = "month", length.out = 12)
  )

  # Capture the error message
  error_message <- tryCatch(
    climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    error = function(e) e$message
  )

  # Check error message
  checkTrue(
    grepl("exactly one value per month", error_message),
    "Error message is not informative for multiple values per month."
  )
}

climdex.pcic.test.SingleMonthlyScalar.raw.different.calendars <- function() {
  set.seed(123)

  # Dates with a "noleap" calendar
  data <- runif(12, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "noleap"), by = "month", length.out = 12)
  result <- try(
    scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
  checkEquals(scalar_obj@dates[!is.na(scalar_obj@data)], dates)
}

climdex.pcic.test.SingleMonthlyScalar.raw.timezones <- function() {
  set.seed(123)

  # Create data with dates including time zones
  data <- runif(12, 0, 20)
  dates <- as.PCICt(seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "month", length.out = 12), cal = "gregorian")
  result <- try(
    scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
  checkEquals(scalar_obj@dates[!is.na(scalar_obj@data)], dates)
}

climdex.pcic.test.SingleMonthlyScalar.raw.irregular.intervals <- function() {
  set.seed(123)

  # Missing months
  data <- runif(10, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "2 months", length.out = 10)
  result <- try(
    scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
    silent = TRUE
  )

  checkTrue(
    !inherits(result, "try-error"),
    "Function raised an error despite valid monthly data."
  )
  checkEquals(scalar_obj@dates[!is.na(scalar_obj@data)], dates)
}

climdex.pcic.test.SingleMonthlyScalar.raw.negative.values <- function() {
  # Data with negative values
  data <- runif(36, -50, 0)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 36)

  scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian")
  scalar_max <- unname(compute.stat.scalar(scalar_obj, "max", "annual", FALSE))
  scalar_min <- unname(compute.stat.scalar(scalar_obj, "min", "annual", FALSE))
  scalar_sum <- unname(compute.stat.scalar(scalar_obj, "sum", "annual", FALSE))
  checkEquals(max(scalar_max), max(data), "annual max stat for single monthly scalar was not equal to max of input data")
  checkEquals(min(scalar_min), min(data), "annual min stat for single monthly scalar was not equal to min of input data")
  checkEquals(sum(scalar_sum), sum(data), "annual sum stat for single monthly scalar was not equal to sum of input data")
}

climdex.pcic.test.SingleMonthlyVector.raw.cartesian <- function() {
  set.seed(123)

  # Data in cartesian format (x,y components)
  x_comp <- runif(12, -10, 10)
  y_comp <- runif(12, -10, 10)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)

  vector_obj <- climdexSingleMonthlyVector.raw(
    primary = x_comp,
    secondary = y_comp,
    dates = dates,
    format = "cartesian",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  # Verify components
  checkEquals(length(vector_obj@primary), length(vector_obj@dates), "Cartesian vector primary length does not match dates.")
  checkEquals(vector_obj@primary[!is.na(vector_obj@primary)], x_comp, "X component not stored correctly")
  checkEquals(vector_obj@secondary[!is.na(vector_obj@secondary)], y_comp, "Y component not stored correctly")
}

climdex.pcic.test.SingleMonthlyScalar.large.dataset <- function() {
  set.seed(123)

  # Create 100 years of monthly data
  n_months <- 100 * 12
  data <- runif(n_months, 0, 20)
  dates <- seq(as.PCICt("1920-01-01", cal = "gregorian"), by = "month", length.out = n_months)

  scalar_obj <- climdexSingleMonthlyScalar.raw(
    data = data,
    dates = dates,
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )

  checkEquals(length(scalar_obj@data[!is.na(scalar_obj@data)]), n_months, "Large dataset not handled correctly")
}
