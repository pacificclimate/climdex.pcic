library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"

# Helper Functions Tests
climdex.pcic.test.convert_cartesian_to_polar <- function() {
  u <- 1
  v <- 1
  result <- convert_cartesian_to_polar(u, v)
  
  checkEqualsNumeric(result$speed, sqrt(2))
  checkEqualsNumeric(result$direction, 45)
  
  checkEqualsNumeric(convert_cartesian_to_polar(0, 1)$direction, 90)
  checkEqualsNumeric(convert_cartesian_to_polar(-1, 0)$direction, 180)
  checkEqualsNumeric(convert_cartesian_to_polar(0, -1)$direction, 270)
  
}

climdex.pcic.test.convert_polar_to_cartesian <- function() {
  speed <- sqrt(2)
  direction <- 45
  result <- convert_polar_to_cartesian(speed, direction)
  
  checkEqualsNumeric(result$u, 1, tolerance = 1e-6)
  checkEqualsNumeric(result$v, 1, tolerance = 1e-6)
  
  checkEqualsNumeric(convert_polar_to_cartesian(0, 45)$u, 0)
  checkEqualsNumeric(convert_polar_to_cartesian(0, 45)$v, 0)
}

climdex.pcic.test.convert_degrees_to_cardinal <- function() {
  # Test all 16 cardinal directions and boundary cases
  degrees <- c(0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 
               202.5, 225, 247.5, 270, 292.5, 315, 337.5, 348.75)
  expected <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 
                'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  result <- convert_degrees_to_cardinal(degrees)
  
  checkEquals(result, expected)
  
  # Boundary cases
  checkEquals(convert_degrees_to_cardinal(c(22.5, 67.5, -45)), c('NNE', 'ENE', 'NW'))
  
  # 360-0 boundary (e.g., 360, should return 'N')
  checkEquals(convert_degrees_to_cardinal(360), 'N')
}


climdex.pcic.test.convert_cardinal_to_degrees <- function() {
  directions <- c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')
  expected <- c(0, 45, 90, 135, 180, 225, 270, 315)
  result <- convert_cardinal_to_degrees(directions)
  
  checkEqualsNumeric(as.numeric(result), expected)
  # Invalid direction should trigger an error
  checkException(convert_cardinal_to_degrees(c('XYZ')), "Invalid cardinal direction provided.")
  
  # Mixed-case input
  mixed_case_result <- convert_cardinal_to_degrees(c('n', 'Ne', 'e', 'SE', 'S', 'sW', 'W', 'nw'))
  checkEqualsNumeric(as.numeric(mixed_case_result), expected)
  
}

climdex.pcic.test.filter_by_direction_range <- function() {
  primary_data <- 1:10
  degrees <- seq(0, 360, length.out = 10)
  date_factors <- rep(1, 10)
  direction.range <- c(90, 270)
  
  filtered <- filter_by_direction_range(primary_data, degrees, date_factors, direction.range)
  
  expected_data <- primary_data[degrees >= 90 & degrees <= 270]
  expected_degrees <- degrees[degrees >= 90 & degrees <= 270]
  
  checkEquals(filtered$primary_data, expected_data)
  checkEquals(filtered$degrees, expected_degrees)
  
  # Test range that crosses 0 
  direction.range <- c(350, 10)
  filtered <- filter_by_direction_range(primary_data, degrees, date_factors, direction.range)
  checkEquals(filtered$primary_data, primary_data[c(1, 10)])
}

climdex.pcic.test.compute_circular_mean <- function() {
  direction_degrees <- c(350, 10)  # Should average to 0 degrees
  date.factors <- c(1,1)
  format <- "degrees"
  result <- compute_circular_mean(direction_degrees, date.factors, format)
  expected <- c(0)
  
  checkEqualsNumeric(result, expected, tolerance = 1e-6)
  
  checkEqualsNumeric(compute_circular_mean(c(90, 90), c(1, 1), "degrees"), 90)
  checkException(compute_circular_mean(c(), c(), "degrees"), 
                 msg = "direction_degrees cannot be empty or NULL.")
}

climdex.pcic.test.compute_circular_sd <- function() {
  direction_degrees <- c(350, 10)
  date.factors <- c(1)
  result <- compute_circular_sd(direction_degrees, date.factors)
  expected <- 10.0261 # Approx
  checkEqualsNumeric(result, expected, tolerance = 1e-4)
  
  # Check for empty inputs (expect an error)
  checkException(compute_circular_sd(c(), c()), 
                 msg = "direction_degrees cannot be empty or NULL.")
}

# # compute.stat.scalar Tests

climdex.pcic.test.compute_stat_scalar_mean <- function() {
  # Test data and dates
  data <- seq(0, 30)
  dates <- seq(as.PCICt("2020-01-01", cal="gregorian"), 
               by = "day", length.out = 31)
  
  # Use the climdexGenericScalar.raw function to create scalar_obj
  scalar_obj <- climdexGenericScalar.raw(
    data = data,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 3, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  # Compute the mean for the 'monthly' frequency
  result <- compute.stat.scalar(scalar_obj, stat = "mean", freq = "monthly", include.exact.dates = FALSE)
  
  expected_mean <- mean(data)
  
  # Compare the computed result for January with the expected value
  checkEqualsNumeric(as.numeric(result[1]), expected_mean)
}


climdex.pcic.test.compute_stat_scalar_max <- function() {
  # Test data and dates
  data <- c(1, 3, 5, 2, 4, 6) 
  dates <- seq(as.PCICt("2020-01-01", cal="gregorian"), by = "day", length.out = 6)
  
  # Use the climdexGenericScalar.raw function to create scalar_obj with higher max.missing.days for monthly
  scalar_obj <- climdexGenericScalar.raw(
    data = data,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 28, seasonal = 6),  # Allow up to 28 missing days in a month
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Compute the maximum value for the 'monthly' frequency
  result <- compute.stat.scalar(scalar_obj, stat = "max", freq = "monthly", include.exact.dates = FALSE)
  
  expected_max <- max(data)
  
  checkEqualsNumeric(as.numeric(result[1]), expected_max)
}
