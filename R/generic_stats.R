library(circular)

#' @title compute.gen.stat
#'
#' @description
#' Computes a specified statistic (e.g., max, min, mean, sum, sd, var) for a generic climate object
#' (either scalar or vector) based on the selected frequency (monthly, annual, or seasonal).
#' The function handles NA masks and exact dates for the max/min statistics.
#'
#' @param gen.var A generic climate object that contains the data, date factors, and NA masks.
#' @param stat The statistic to compute. Must be one of `"max"`, `"min"`, `"mean"`, `"sum"`, `"sd"`, or `"var"`.
#' @param data The numeric data on which to compute the statistic.
#' @param freq The frequency for which to compute the statistic. Options are `"monthly"`, `"annual"`, or `"seasonal"`.
#' @param include.exact.dates Logical. If `TRUE`, computes and returns the exact dates for max/min statistics.
#' This is not applicable for other statistics like mean, sum, sd, or var, and will automatically switch to `FALSE`.
#'
#' @return A numeric vector of the computed statistic for each date factor, multiplied by the corresponding NA mask.
#' If `include.exact.dates = TRUE`, it returns a data frame that includes exact dates for max/min statistics.
#'
#' @details
#' The function checks the validity of the provided data and ensures that NA masks are applied to the results.
#' If exact dates are requested for statistics that do not support it (e.g., mean, sum, sd, var), the function will print a message
#' and proceed without exact dates.
#'
#' @seealso \code{\link{compute.stat.scalar}}, \code{\link{compute.stat.vector}}
#'
#' @note
#' This function is internal and not intended to be called directly by users. It serves as a shared utility
#' for computing statistics of generic scalar and vector climate data.
#'
#' @examples
#' # Assuming `scalar_obj` is a valid climdexGenericScalar object:
#' \dontrun{compute.gen.stat(scalar_obj, "max", scalar_obj@data, "monthly", FALSE)}
#'
#'#' @export
compute.gen.stat <- function(gen.var, stat, data, freq = c("monthly", "annual", "seasonal"), include.exact.dates) {
  stopifnot(!is.null(data))
  freq <- match.arg(freq)
  exact_date_stats <- c("max", "min")
  
  if (include.exact.dates && !(stat %in% exact_date_stats)) {
    message(paste("Warning: Exact dates are not applicable for the", stat, "statistic. Proceeding without exact dates."))
    include.exact.dates <- FALSE
  }
  
  date.factors <- gen.var@date.factors[[freq]]
  mask <- gen.var@namasks[[freq]][[1]]
  cal <- attr(gen.var@dates, "cal")

  result <- suppressWarnings(tapply.fast(data, date.factors, stat, na.rm = TRUE))
 
  if (include.exact.dates) {
    return(exact.date(stat, data, date.factors, freq, cal, mask))
  }

  return(result * mask)
}

#' @title compute.stat.scalar
#'
#' @description
#' Computes a specified statistic (e.g., max, min, mean, sum, sd, var) for scalar climate data
#' based on a given frequency (monthly, annual, seasonal). It can also compute exact dates for
#' specific statistics when requested.
#'
#' @param scalar_obj A `ClimdexGenericScalar` object containing the scalar climate data.
#' @param stat The statistic to compute. Options are `"max"`, `"min"`, `"mean"`, `"sum"`, `"sd"`, or `"var"`.
#' @param freq The frequency of the statistic. Options are `"monthly"`, `"annual"`, or `"seasonal"`.
#' @param include.exact.dates Logical. If `TRUE`, returns the exact dates of the max/min values. Not applicable for `mean`, `sd`, or `var`.
#'
#' @return The computed statistic for the given frequency and statistic type.
#'
#' @seealso \code{\link{compute.gen.stat}}
#'
#' @examples
#' # Example usage for scalar data:
#' \dontrun{compute.stat.scalar(scalar_obj, "max", "monthly", TRUE)}
#'
#' @export
compute.stat.scalar <- function(scalar_obj, 
                                stat = c("max", "min", "mean", "sum", "sd", "var"),
                                freq = c("monthly", "annual", "seasonal"), include.exact.dates = FALSE
) {
  stat <- match.arg(stat)
  stopifnot(!is.null(scalar_obj@data))  # Ensure the data key exists
  return(compute.gen.stat(scalar_obj, stat, scalar_obj@data, freq, include.exact.dates))
}

#' @title Convert Cartesian to Polar Coordinates
#' @description Converts Cartesian coordinates (u, v) to polar coordinates (speed, direction).
#' @param u A numeric value representing the x-component (u) of the Cartesian coordinate.
#' @param v A numeric value representing the y-component (v) of the Cartesian coordinate.
#' @return A list with two elements:
#'   \item{speed}{The magnitude of the vector (i.e., the speed).}
#'   \item{direction}{The direction of the vector in degrees, normalized to the range [0, 360).}
#' @export
convert_cartesian_to_polar <- function(u, v) {
  if (!is.numeric(u) || !is.numeric(v)) {
    stop("u and v must be numeric.")
  }
  speed <- sqrt(u^2 + v^2)
  
  if (speed == 0) {
    direction <- NA  # No meaningful direction when speed is 0
  } else {
    direction <- (atan2(v, u) * 180 / pi + 360) %% 360
  }
  
  direction <- (atan2(v, u) * 180 / pi + 360) %% 360
  return(list(speed = speed, direction = direction))
}

#' @title Convert Polar to Cartesian Coordinates
#' @description Converts polar coordinates (speed, direction) to Cartesian coordinates (u, v).
#' @param speed A numeric value representing the magnitude of the vector.
#' @param direction A numeric value representing the direction of the vector in degrees.
#' @return A list with two elements:
#'   \item{u}{The x-component (u) of the Cartesian coordinate.}
#'   \item{v}{The y-component (v) of the Cartesian coordinate.}
#' @export
convert_polar_to_cartesian <- function(speed, direction) {
  if (!is.numeric(speed) || !is.numeric(direction)) {
    stop("speed and direction must be numeric.")
  }
  radians <- direction * pi / 180
  u <- speed * cos(radians)
  v <- speed * sin(radians)
  return(list(u = u, v = v))
}

#' @title Convert Degrees to Cardinal Direction
#' @description Converts a numeric degree value into the corresponding cardinal direction (e.g., N, NE, E).
#' @param degrees A numeric vector of degrees (0-360) representing the direction.
#' @return A character vector of cardinal directions corresponding to the degree values.
#' @export
convert_degrees_to_cardinal <- function(degrees) {
  if (!is.numeric(degrees)) {
    stop("degrees must be a numeric vector.")
  }
  directions <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 
                  'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  index <- round(degrees / 22.5) %% 16 + 1
  return(directions[index])
}

#' @title Convert Cardinal Direction to Degrees
#' @description Converts a cardinal direction (e.g., N, NE, E) into the corresponding degree value.
#' @param direction A character vector of cardinal directions.
#' @return A numeric vector of degrees corresponding to the cardinal directions.
#' @export
convert_cardinal_to_degrees <- function(direction) {
  cardinal_map <- c(N = 0, NNE = 22.5, NE = 45, ENE = 67.5, E = 90,
                    ESE = 112.5, SE = 135, SSE = 157.5, S = 180,
                    SSW = 202.5, SW = 225, WSW = 247.5, W = 270,
                    WNW = 292.5, NW = 315, NNW = 337.5)
  
  direction <- toupper(direction)
  if (any(!direction %in% names(cardinal_map))) {
    stop("Invalid cardinal direction provided.")
  }
  return(cardinal_map[direction])
}

#' @title Filter Data by Direction Range
#' @description Filters primary data, degrees, and date factors based on a specified range of directions.
#' @param primary_data A numeric vector of the primary data to be filtered.
#' @param degrees A numeric vector of direction degrees.
#' @param date_factors A vector of date factors corresponding to the data.
#' @param direction.range A numeric vector of length 2 specifying the minimum and maximum degrees for filtering.
#' @return A list containing:
#'   \item{primary_data}{The filtered primary data.}
#'   \item{degrees}{The filtered degrees.}
#'   \item{date_factors}{The filtered date factors.}
#' @export
filter_by_direction_range <- function(primary_data, degrees, date_factors, direction.range) {
  if (!is.numeric(direction.range) || length(direction.range) != 2) {
    stop("direction.range must be a numeric vector of length 2 specifying min and max degrees.")
  }
  
  # Normalize degrees to [0, 360)
  degrees_normalized <- (degrees + 360) %% 360
  min_dir <- direction.range[1] %% 360
  max_dir <- direction.range[2] %% 360
  
  # Handle ranges that cross the 0-degree line
  if (min_dir > max_dir) {
    within_range <- degrees_normalized >= min_dir | degrees_normalized <= max_dir
  } else {
    within_range <- degrees_normalized >= min_dir & degrees_normalized <= max_dir
  }
  
  # Apply the filter to data and date factors
  primary_data_filtered <- primary_data[within_range]
  degrees_filtered <- degrees[within_range]
  date_factors_filtered <- date_factors[within_range]
  
  return(list(
    primary_data = primary_data_filtered,
    degrees = degrees_filtered,
    date_factors = date_factors_filtered
  ))
}

#' @title Compute Circular Mean
#' @description Computes the circular mean of directional data based on date factors.
#' @param direction_degrees A numeric vector of direction values in degrees.
#' @param date.factors A vector of date factors used to group the data.
#' @param format The format of the result, either "degrees" (default) or "cardinal".
#' @return A numeric vector of circular means for each date factor group, or a character vector of cardinal directions if `format` is set to "cardinal".
#' @importFrom circular circular mean.circular sd.circular
#' @export
compute_circular_mean <- function(direction_degrees, date.factors, format) {
  
  # Check for empty or NULL inputs
  if (is.null(direction_degrees) || length(direction_degrees) == 0) {
    stop("direction_degrees cannot be empty or NULL.")
  }
  if (is.null(date.factors) || length(date.factors) == 0) {
    stop("date.factors cannot be empty or NULL.")
  }
  
  # Check if inputs are numeric
  if (!is.numeric(direction_degrees)) {
    stop("direction_degrees must be a numeric vector.")
  }
  
  # Check if format is valid
  if (!format %in% c("degrees", "cardinal")) {
    stop("format must be either 'degrees' or 'cardinal'.")
  }
  
  # Convert directions to 'circular' objects
  directions_circular <- circular::circular(direction_degrees, units = "degrees", modulo = "2pi")
  
  # Compute circular mean
  circular_mean <- tapply(directions_circular, date.factors, circular::mean.circular, na.rm = TRUE)
  
  # Convert back to degrees
  circular_mean_degrees <- as.numeric(circular_mean)
  circular_mean_degrees <- (circular_mean_degrees + 360) %% 360  # Normalize to [0, 360)
  
  # Convert to cardinal if format is 'cardinal'
  if (format == "cardinal") {
    circular_mean_result <- sapply(circular_mean_degrees, convert_degrees_to_cardinal)
  } else {
    circular_mean_result <- circular_mean_degrees
  }
  
  return(circular_mean_result)
}

#' @title Compute Circular Standard Deviation
#' @description Computes the circular standard deviation of directional data based on date factors.
#' @param direction_degrees A numeric vector of direction values in degrees.
#' @param date.factors A vector of date factors used to group the data.
#' @return A numeric vector of circular standard deviations for each date factor group, in degrees.
#' @export
compute_circular_sd <- function(direction_degrees, date.factors) {

  # Check for empty or NULL inputs
  if (is.null(direction_degrees) || length(direction_degrees) == 0) {
    stop("direction_degrees cannot be empty or NULL.")
  }
  if (is.null(date.factors) || length(date.factors) == 0) {
    stop("date.factors cannot be empty or NULL.")
  }
  
  # Check if inputs are numeric
  if (!is.numeric(direction_degrees)) {
    stop("direction_degrees must be a numeric vector.")
  }

  # Convert directions to 'circular' objects
  directions_circular <- circular::circular(direction_degrees, units = "degrees", modulo = "2pi")
  
  # Compute circular standard deviation
  circular_sd <- tapply(directions_circular, date.factors, circular::sd.circular, na.rm = TRUE)
  circular_sd_degrees <- as.numeric(circular_sd) * (180 / pi)  # Convert from radians to degrees
  
  return(circular_sd_degrees)
}

#' @title compute.stat.vector
#'
#' @description
#' Computes a specified statistic (e.g., max, min, mean, sd, circular_mean, etc.) for vector climate data.
#' The data can be in polar, cartesian, or cardinal format, and the function handles both magnitude and direction.
#' Supports additional vector-specific statistics like circular mean and circular standard deviation.
#'
#' @param climate_obj A `ClimdexGenericVector` object containing vector climate data (e.g., wind speed and direction).
#' @param stat The statistic to compute. Must be one of `"max"`, `"min"`, `"mean"`, `"sum"`, `"circular_mean"`, `"sd"`, `"var"`, or `"circular_sd"`.
#' @param freq The frequency for which to compute the statistic. Options are `"monthly"`, `"annual"`, or `"seasonal"`.
#' @param format The format of the vector data. Must be one of `"polar"`, `"cartesian"`, or `"cardinal"`.
#' This determines how the primary and secondary components (e.g., magnitude and direction) are interpreted.
#' @param include.exact.dates Logical. If `TRUE`, computes and returns the exact dates for max/min statistics. Not applicable for other statistics.
#' @param direction.range A numeric vector of length 2 specifying the minimum and maximum degrees for filtering based on direction.
#' Only data with directions within this range will be included in the calculations. If `NULL`, no filtering is applied.
#'
#' @return A list containing the computed statistic for magnitude and, if applicable, the computed statistic for direction.
#' For circular statistics (e.g., `circular_mean`, `circular_sd`), the result is returned for directions in degrees.
#' If `include.exact.dates = TRUE`, the function returns exact dates for the max/min statistics.
#'
#' @details
#' The function computes the specified statistic for both the magnitude (primary component) and direction (secondary component) of the vector data.
#' It supports additional statistics like `circular_mean` and `circular_sd` for directional data.
#' If `include.exact.dates = TRUE`, exact dates are returned for max/min statistics.
#' The function can also filter the data based on a specified range of directions (using `direction.range`).
#'
#' @seealso \code{\link{compute.stat.scalar}}, \code{\link{compute.gen.stat}}
#'
#' @note
#' This function is designed for vector climate data, where the data includes both a magnitude and direction component.
#' For scalar data, use \code{\link{compute.stat.scalar}} instead.
#'
#' @examples
#' \dontrun{
#' # Assuming `vector_obj` is a valid ClimdexGenericVector object:
#' compute.stat.vector(vector_obj, "circular_mean", "monthly", format = "polar")
#'}
#'
compute.stat.vector <- function(
    climate_obj,
    stat = c("max", "min", "mean", "sum", "circular_mean", "sd", "var", "circular_sd"),
    freq = c("monthly", "annual", "seasonal"),
    format = c("polar", "cartesian", "cardinal"),
    include.exact.dates = FALSE,
    direction.range = NULL
) {
  stat <- match.arg(stat)
  freq <- match.arg(freq)
  format <- match.arg(format)
  
  date.factors <- climate_obj@date.factors[[freq]]
  
  # Convert all data to polar coordinates (magnitude and direction in degrees)
  switch(
    format,
    "cartesian" = {
      polar_data <- convert_cartesian_to_polar(climate_obj@primary, climate_obj@secondary)
      magnitude <- polar_data$speed
      direction_degrees <- polar_data$direction
    },
    "cardinal" = {
      direction_degrees <- sapply(climate_obj@secondary, convert_cardinal_to_degrees)
      magnitude <- climate_obj@primary
    },
    "polar" = {
      magnitude <- climate_obj@primary
      direction_degrees <- as.numeric(climate_obj@secondary)
    }
  )
  
  # Filter data based on direction range if provided
  if (!is.null(direction.range)) {
    filtered <- filter_by_direction_range(magnitude, direction_degrees, date.factors, direction.range)
    magnitude <- filtered$primary_data
    direction_degrees <- filtered$degrees
    date.factors <- filtered$date_factors
  }
  
  result <- switch(
    stat,
    "circular_mean" = {
      direction_result <- compute_circular_mean(direction_degrees, date.factors, format)
      list(direction = direction_result)
    },
    "circular_sd" = {
      circular_sd_degrees <- compute_circular_sd(direction_degrees, date.factors)
      list(circular_sd = circular_sd_degrees)
    },
    {
      # For other statistics, compute on magnitude
      # Create a local copy of climate_obj to avoid modifying the original object
      local_climate_obj <- climate_obj
      
      # Define a unique key for magnitude data
      magnitude_key <- "magnitude_temp"
      
      # Add magnitude data to local_climate_obj@data
      local_climate_obj@data[[magnitude_key]] <- magnitude
      
      # Use existing masks for primary and secondary data
      primary <- climate_obj@primary
      secondary <- climate_obj@secondary
      mask_primary <- climate_obj@namasks[[freq]][[primary]]
      mask_secondary <- climate_obj@namasks[[freq]][[secondary]]
      
      if (!is.null(mask_primary) && !is.null(mask_secondary)) {
        magnitude_mask <- mask_primary * mask_secondary
      } else if (!is.null(mask_primary)) {
        magnitude_mask <- mask_primary
      } else if (!is.null(mask_secondary)) {
        magnitude_mask <- mask_secondary
      } else {
        # Default mask
        magnitude_mask <- rep(1, length(unique(date.factors)))
        names(magnitude_mask) <- unique(date.factors)
      }
      
      # Add the mask to local_climate_obj@namasks
      local_climate_obj@namasks[[freq]][[magnitude_key]] <- magnitude_mask
      
      magnitude_stat <- compute.gen.stat(
        local_climate_obj,
        stat,
        magnitude_key,
        freq = freq,
        include.exact.dates = include.exact.dates
      )
      
      list(magnitude = magnitude_stat)
    }
  )
  
  return(result)
}

