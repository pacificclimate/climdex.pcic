#' @title climdexGenericVector.raw
#' 
#' @description
#' Creates a `ClimdexGenericVector` object from raw vector climate data, including
#' both a primary (e.g., magnitude) and secondary (e.g., direction) component.
#' 
#' @details
#' This function processes vector climate data and creates a `ClimdexGenericVector`
#' object. The function generates NA masks based on the provided `max.missing.days`
#' and validates the `primary` and `secondary` components based on the specified
#' format (`polar`, `cartesian`, or `cardinal`).
#' 
#' @param primary A numeric vector representing the primary data (e.g., wind speed).
#' @param secondary A numeric or character vector representing the secondary data (e.g., wind direction).
#' @param dates A `PCICt` vector corresponding to the data dates.
#' @param format A string specifying the format of the vector data ("polar", "cartesian", "cardinal").
#' @param max.missing.days A named vector indicating the maximum allowed missing days for `annual`, `monthly`, and `seasonal` time periods.
#' @param northern.hemisphere Whether this point is in the northern hemisphere.
#' @param calendar String representing the calendar type, e.g., "gregorian".
#' 
#' @return A `ClimdexGenericVector` object containing the processed vector data.
#' 
#' @seealso \code{\link{ClimdexGenericScalar.raw}}, \code{\link{ClimdexGenericVector.csv}}
#' 
#' @examples
#' primary <- c(5.5, 6.2, 4.8)
#' secondary <- c("N", "NE", "E")
#' dates <- as.PCICt(c("2000-01-01", "2000-01-02", "2000-01-03"), format = "%Y-%m-%d", cal = "gregorian")
#' climdexGenericVector.raw(primary, secondary, dates, format = "cardinal")
#' 
#' @export

climdexGenericVector.raw <- function(
    primary,
    secondary,
    dates,
    format = "polar",
    max.missing.days = c(annual = 15, monthly = 3, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
) {

  check.generic.argument.validity(list(primary, secondary), dates, calendar)
  
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
  
  date.info <- date_info(dates)
  jdays = date.info$jdays
  date.series = date.info$date.series
  date.factors = date.info$date.factors

  filled.list <- generate_filled_list(list(primary, secondary), dates, date.series)
  filled.primary <- filled.list[[1]]
  filled.secondary <- filled.list[[2]]

  namasks <- generate_namasks(filled.list, date.factors, max.missing.days)

  obj <- new("ClimdexGenericVector",
    primary = filled.primary,
    secondary = filled.secondary,
    dates = date.series,
    format = format,
    date.factors = date.factors,
    namasks = namasks,
    max.missing.days = max.missing.days,
    northern.hemisphere = northern.hemisphere)
  
  return(obj)
}

#' @title climdexGenericVector.csv
#'
#' @description
#' Reads vector climate data (e.g., wind speed and direction) from a CSV file and creates a `ClimdexGenericVector` object.
#'
#' @details
#' This function reads vector climate data, consisting of primary and secondary components (e.g., magnitude and direction),
#' from a CSV file. It validates the data and constructs a `ClimdexGenericVector` object that can be used to compute
#' climate indices or for further analysis. The format of the vector data (e.g., `polar`, `cartesian`, or `cardinal`)
#' must be specified.
#'
#' @param file The file path to the CSV containing the vector climate data.
#' @param primary.column The name of the column containing the primary data (e.g., magnitude) in the CSV file.
#' @param secondary.column The name of the column containing the secondary data (e.g., direction) in the CSV file.
#' @param date.columns A vector of column names or indices corresponding to the date fields in the CSV file.
#' @param date.format A string representing the format of the date fields.
#' @param name A name for the vector data (for reference purposes).
#' @param format A string specifying the format of the vector data. Must be one of `"polar"`, `"cartesian"`, or `"cardinal"`.
#' @param na.strings A character vector of strings to interpret as `NA`.
#' @param max.missing.days A named vector specifying the maximum number of missing days allowed for annual, monthly, and seasonal periods.
#' @param northern.hemisphere Whether this point is in the northern hemisphere.
#' @param calendar A string representing the calendar type (e.g., "gregorian").
#'
#' @return A `ClimdexGenericVector` object containing the processed vector climate data.
#'
#' @seealso \code{\link{ClimdexGenericVector.raw}}, \code{\link{ClimdexGenericScalar.csv}}
#'
#' @examples
#' # Example usage:
#' climdexGenericVector.csv("data.csv", primary.column = "wind_speed", secondary.column = "wind_direction",
#'                          date.columns = c("year", "month", "day"), date.format = "%Y %m %d",
#'                          format = "cardinal", calendar = "gregorian")
#'
#' @export

climdexGenericVector.csv <- function(
    file,
    primary.column,
    secondary.column,
    date.columns,
    date.format,
    name,
    format = "polar",
    na.strings = NULL,
    max.missing.days = c(annual = 15, monthly = 3, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
) {

  GV.csv <- read_csv_data(file, data.columns = c(primary.column, secondary.column), date.columns, date.format, na.strings, calendar)
  
  primary_values <- GV.csv$data[[1]]
  secondary_values <- GV.csv$data[[2]]
  dates <- GV.csv$dates
  
  obj <- climdexGenericVector.raw(
    primary = primary_values,
    secondary = secondary_values,
    dates = dates,
    format = format,
    max.missing.days = max.missing.days,
    northern.hemisphere = northern.hemisphere,
    calendar = calendar
  )
  
  return(obj)
}