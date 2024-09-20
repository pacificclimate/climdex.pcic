#' Class definition for ClimdexGenericVariable
#' @title ClimdexGenericVariable
#' 
#' @description
#' This class represents the base for handling both scalar and vector generic climate data.
#' It is used internally to construct single-variable, climdexInput-like objects for caclulating basic climate indices.
#' 
#' @section Slots:
#' \describe{
#'   \item{dates}{A `PCICt` type vector of dates.}
#'   \item{date.factors}{Factors used for creation of annual, seasonal and monthly indices.}
#'   \item{jdays}{Julian days for the date sequence.}
#'   \item{namasks}{Data quality masks for annual, seasonal and monthly data.}
#'   \item{northern.hemisphere}{Logical indicating whether data is from the northern hemisphere.}
#'   \item{max.missing.days}{A named vector specifying the maximum number of missing days of data for annual, seasonal
#' and monthly data.}
#' }
#' 
#' @keywords internal

setClass(
  "ClimdexGenericVariable",
  slots = list(
    dates = "PCICt",
    date.factors = "list",
    jdays="numeric",
    namasks = "list",
    northern.hemisphere = "logical",
    max.missing.days = "numeric"
  )
)

## Class definition declaration
#' @title ClimdexGenericScalar
#' 
#' @description
#' The `ClimdexGenericScalar` class contains scalar climate data (e.g., humidity, snow-depth
#' ) and other data required for days for calculating basic climate indices.
#'
#' @section Slots:
#' \describe{
#'   \item{data}{A numeric vector of climate data.}
#'   \item{dates}{A `PCICt` type vector of dates.}
#'   \item{date.factors}{Factors used for creation of annual, seasonal and monthly indices.}
#'   \item{jdays}{Julian days for the date sequence.}
#'   \item{namasks}{Data quality masks for annual, seasonal and monthly data.}
#'   \item{northern.hemisphere}{Logical indicating whether data is from the northern hemisphere.}
#'   \item{max.missing.days}{A named vector specifying the maximum number of missing days of data for annual, seasonal
#' and monthly data.}
#' }
#' 
#' @name ClimdexGenericScalar
#' @aliases ClimdexGenericScalar-class
#' @exportClass ClimdexGenericScalar

setClass(
  "ClimdexGenericScalar",
  contains = "ClimdexGenericVariable",
  slots = list(
    data = "numeric"
  )
)

## Class definition declaration
#' @title ClimdexGenericVector
#' 
#' @description
#' The `ClimdexGenericVector` class contains vector climate data (e.g., wind speed and direction), consisting of
#' both primary and secondary components, and other data required for days for calculating basic climate indices.
#'
#' @details
#' The `primary` slot contains the main numeric climate data, such as wind speed or another magnitude-based variable.
#' The `secondary` slot provides complementary information, such as direction, and its type and meaning depend on the 
#' `format`:
#' \itemize{
#'   \item \strong{polar}: When the format is `"polar"`, the `primary` slot represents magnitude (e.g., wind speed),
#'   and the `secondary` slot must be numeric and represents direction in degrees (e.g., wind direction, where 0° is North).
#'   \item \strong{cartesian}: In `"cartesian"` format, `primary` and `secondary` both represent numeric values, typically the
#'   x and y components of a vector (e.g., velocity components in the x and y directions).
#'   \item \strong{cardinal}: When the format is `"cardinal"`, the `primary` slot still represents a numeric value (e.g., wind speed),
#'   but the `secondary` slot is a character vector representing cardinal directions (e.g., "N", "NE", "E", etc.).
#' }
#' This class is designed to handle both magnitude and direction data in various forms, allowing flexibility in representing
#' vector climate data depending on the source.
#'
#' @section Slots:
#' \describe{
#'   \item{primary}{A numeric vector representing the primary climate data (e.g., wind speed).}
#'   \item{secondary}{A vector representing the secondary climate data (e.g., wind direction).}
#'   \item{format}{A string indicating the format of the data: "polar", "cartesian", or "cardinal".}
#'   \item{dates}{A `PCICt` type vector of dates.}
#'   \item{date.factors}{Factors used for creation of annual, seasonal and monthly indices.}
#'   \item{jdays}{Julian days for the date sequence.}
#'   \item{namasks}{Data quality masks for annual, seasonal and monthly data.}
#'   \item{northern.hemisphere}{Logical indicating whether data is from the northern hemisphere.}
#'   \item{max.missing.days}{A named vector specifying the maximum number of missing days of data for annual, seasonal
#'    and monthly data.}
#' }
#' 
#' @name ClimdexGenericVector
#' @aliases ClimdexGenericVector-class
#' @exportClass ClimdexGenericVector

setClass(
  "ClimdexGenericVector",
  contains = "ClimdexGenericVariable",
  slots = list(
    primary = "numeric",
    secondary = "ANY",  # Could be numeric or character, depending on format.
    format = "character"  # 'polar', 'cartesian', or 'cardinal'
  )
)
