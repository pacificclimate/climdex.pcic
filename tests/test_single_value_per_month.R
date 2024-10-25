library(climdex.pcic)
library(RUnit)

climdex.pcic.test.single.monthly.scalar.raw.and.csv.construction <- function() {
  set.seed(123)
  
  # Raw data construction using 1st of each month
  scalar_data <- runif(12, 0, 20)  # One value per month
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)
  scalar_obj_raw <- climdexSingleMonthlyScalar.raw(
    data = scalar_data,
    dates = dates,
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  checkEquals(length(scalar_obj_raw@data), length(scalar_obj_raw@dates), "Raw scalar construction: data length does not match dates.")
  
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
  
  checkEquals(length(scalar_obj_csv@data), length(scalar_obj_csv@dates), "CSV scalar construction: data length does not match dates.")
  checkEquals(scalar_obj_raw@dates, scalar_obj_csv@dates, "Date mismatch between raw and CSV scalar objects.")
  checkTrue(all.equal(scalar_obj_csv, scalar_obj_raw), msg = "Scalar_obj built from CSV is not identical to raw")
}

climdex.pcic.test.single.monthly.vector.raw.and.csv.construction <- function() {
  set.seed(123)
  
  # Raw data construction using 1st of each month
  primary_data <- runif(12, 0, 20)  # One value per month
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
  
  checkEquals(length(vector_obj_raw@primary), length(vector_obj_raw@dates), "Raw vector construction: primary length does not match dates.")
  
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
  
  checkEquals(length(vector_obj_csv@primary), length(vector_obj_csv@dates), "CSV vector construction: primary length does not match dates.")
  checkEquals(vector_obj_raw@dates, vector_obj_csv@dates, "Date mismatch between raw and CSV vector objects.")
  checkTrue(all.equal(vector_obj_csv, vector_obj_raw), msg = "Vector_obj built from CSV is not identical to raw")
}

# Test for scalar with missing values in data
climdex.pcic.test.SingleMonthlyScalar.raw.missing <- function() {
  set.seed(123)
  
  # Simulate single monthly value data with an NA value
  data <- c(1:5, NA, 7:12)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)
  
  # Create the climdexSingleMonthlyScalar object and expect it to pass without failure
  scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian")
  
  checkTrue(length(scalar_obj@data) == length(scalar_obj@dates), msg = "NA values in single monthly scalar data are not handled correctly.")
}

# Test for vector with missing values in primary data
climdex.pcic.test.SingleMonthlyVector.raw.missing <- function() {
  set.seed(123)
  
  # Simulate single monthly value vector data with an NA value in the primary component
  primary <- c(1:5, NA, 7:12)
  secondary <- runif(12, 0, 360)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "month", length.out = 12)
  
  # Create the climdexSingleMonthlyVector object and expect it to pass without failure
  vector_obj <- climdexSingleMonthlyVector.raw(primary, secondary, dates, format = "polar", northern.hemisphere = TRUE, calendar = "gregorian")
  
  checkTrue(length(vector_obj@primary) == length(vector_obj@dates), msg = "NA values in single monthly vector data are not handled correctly.")
}

# Test with zero-length scalar data
climdex.pcic.test.SingleMonthlyScalar.raw.zero.length <- function() {
  set.seed(123)
  
  # Create an empty vector and dates
  data <- numeric(0)
  dates <- as.PCICt(character(0), cal = "gregorian")
  checkException(climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian"),
                 "Zero-length data not handled correctly.")
}

# Test with out-of-sequence dates and sorted check
climdex.pcic.test.OutOfSequenceDates <- function() {
  set.seed(123)
  
  # Create scalar data with out-of-sequence dates
  data <- 1:12  
  dates <- as.PCICt(c("2020-01-01", "2020-03-01", "2020-05-01", "2020-02-01", "2020-04-01", "2020-06-01",
                      "2020-07-01", "2020-08-01", "2020-09-01", "2020-10-01", "2020-11-01", "2020-12-01"),
                    cal = "gregorian")
  
  # Create the scalar object with out-of-sequence dates
  scalar_obj <- climdexSingleMonthlyScalar.raw(data, dates, northern.hemisphere = TRUE, calendar = "gregorian")
  
  # Extract the non-NA data and corresponding dates from the object
  obj_dates <- scalar_obj@dates[!is.na(scalar_obj@data)]
  obj_data <- scalar_obj@data[!is.na(scalar_obj@data)]
  
  sorted_indices <- order(dates)
  sorted_dates <- dates[sorted_indices]
  sorted_data <- data[sorted_indices]
  
  checkEquals(obj_data, sorted_data, 
              "Object Data is not aligned correctly with sorted dates.")
}
