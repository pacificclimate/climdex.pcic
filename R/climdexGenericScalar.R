#' @title climdexGenericScalar.raw
#' 
#' @description
#' Creates a `ClimdexGenericScalar` object from raw scalar climate data.
#' 
#' @details
#' This function processes scalar climate data (e.g., humidity, snow-depth)
#' and creates a `ClimdexGenericScalar` object. The function fills missing
#' values and applies NA masks based on the `max.missing.days` argument.
#' The `ClimdexGenericScalar` class is used to compute
#' basic climate indices from scalar data.
#' 
#' @param data A numeric vector containing the scalar climate data.
#' @param dates A `PCICt` vector corresponding to the data dates.
#' @param max.missing.days A named vector indicating the maximum allowed missing days for `annual`, `monthly`, and `seasonal` time periods.
#' @param northern.hemisphere Whether this point is in the northern hemisphere.
#' @param calendar String representing the calendar type, e.g., "gregorian".
#' @return A `ClimdexGenericScalar` object containing the processed data.
#' 
#' @seealso \code{\link{climdexGenericVector.raw}}, \code{\link{climdexGenericScalar.csv}}
#' 
#' @examples
#' data <- c(10.5, 12.3, 11.2)
#' dates <- as.PCICt(c("2024-01-01", "2024-01-02", "2024-01-03"),
#'                   format = "%Y-%m-%d", cal = "gregorian")
#' climdexGenericScalar.raw(data, 
#'                          dates,
#'                          max.missing.days = c(annual = 15, monthly = 3, seasonal = 6))
#' 
#' @export

climdexGenericScalar.raw <- function(
  data,
  dates,
  max.missing.days = c(annual = 15, monthly = 3, seasonal = 6),
  northern.hemisphere = TRUE,
  calendar = "gregorian"
) {
  
  check.generic.argument.validity(data,dates,max.missing.days)

  date.info <- date_info(dates)
  jdays = date.info$jdays
  date.series = date.info$date.series
  date.factors = date.info$date.factors

  filled.list <- generate_filled_list(data, dates, date.series)
  names(filled.list) <- "data"
  namasks <- generate_namasks(filled.list, date.factors, max.missing.days)
  obj <- new("climdexGenericScalar",
    data = filled.list[["data"]],
    dates = date.series,
    date.factors = date.factors,
    jdays = jdays,
    namasks = namasks,
    northern.hemisphere = northern.hemisphere,
    max.missing.days = max.missing.days)
  
  return(obj)
}

#' @title climdexGenericScalar.csv
#'
#' @description
#' Reads scalar climate data from a CSV file and creates a `ClimdexGenericScalar` object.
#'
#' @details
#' This function reads scalar climate data (e.g., humidity, snow-depth) from a CSV file, validates the data,
#' and generates a `ClimdexGenericScalar` object.
#'
#' The CSV file should contain the climate data in the specified columns, and the date fields should be provided in separate columns.
#'
#' @param file The file path to the CSV containing the scalar climate data.
#' @param data.column The name of the column containing the scalar data in the CSV file.
#' @param date.columns A vector of column names corresponding to the date fields in the CSV file.
#' @param date.format A string representing the format of the date fields.
#' @param na.strings A character vector of strings to interpret as `NA`.
#' @param northern.hemisphere Logical indicating whether the data is from the northern hemisphere.
#' @param max.missing.days A named vector specifying the maximum number of missing days allowed for annual, monthly, and seasonal periods.
#' @param calendar A string representing the calendar type (e.g., "gregorian").
#'
#' @return A `ClimdexGenericScalar` object containing the processed scalar climate data.
#'
#' @seealso \code{\link{climdexGenericScalar.raw}}, \code{\link{climdexGenericVector.csv}}
#'
#' @examples
#' # Example usage for scalar data:
#'
#' # Simulating CSV data for humidity
#' csv_data <- "
#' year,month,day,humidity
#' 2024,01,01,80
#' 2024,01,02,82
#' 2024,01,03,85
#' "
#'
#' # Write the CSV to a temporary file
#' temp_file <- tempfile(fileext = ".csv")
#' writeLines(csv_data, temp_file)  
#'
#' # Call the climdexGenericScalar.csv function
#' climdexGenericScalar.csv(temp_file, data.column = "humidity",
#'                          date.columns = c("year", "month", "day"),
#'                          date.format = "%Y %m %d", calendar = "gregorian")

#' @export

climdexGenericScalar.csv <- function(
    file,
    data.column,
    date.columns,
    date.format,
    na.strings = NULL,
    northern.hemisphere = TRUE,
    max.missing.days = c(annual = 15, monthly = 3, seasonal = 6),
    calendar = "gregorian"
) {

  GS.csv <- read_csv_data(file, data.column, date.columns, date.format,  na.strings, calendar)
  obj <- climdexGenericScalar.raw(
    data = GS.csv$data[[1]],
    dates = GS.csv$dates,
    northern.hemisphere = northern.hemisphere,
    max.missing.days = max.missing.days,
    calendar = calendar
  )
  
  return(obj)
}
