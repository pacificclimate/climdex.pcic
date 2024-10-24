library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"

# Helper Functions Tests
climdex.pcic.test.convert.cartesian.to.polar <- function() {
  u <- 1
  v <- 1
  result <- convert_cartesian_to_polar(u, v)
  
  checkEqualsNumeric(result$speed, sqrt(2))
  checkEqualsNumeric(result$direction, 45)
  
  checkEqualsNumeric(convert_cartesian_to_polar(0, 1)$direction, 90)
  checkEqualsNumeric(convert_cartesian_to_polar(-1, 0)$direction, 180)
  checkEqualsNumeric(convert_cartesian_to_polar(0, -1)$direction, 270)
  
  result_zero <- convert_cartesian_to_polar(0, 0)
  checkEqualsNumeric(result_zero$speed, 0, "Expected speed to be 0 when both u and v are 0")
  checkTrue(is.na(result_zero$direction), "Expected direction to be NA when both u and v are 0")
}

climdex.pcic.test.convert.polar.to.cartesian <- function() {
  speed <- sqrt(2)
  direction <- 45
  result <- convert_polar_to_cartesian(speed, direction)
  
  checkEqualsNumeric(result$u, 1, tolerance = 1e-6)
  checkEqualsNumeric(result$v, 1, tolerance = 1e-6)
  
  # No speed
  checkEqualsNumeric(convert_polar_to_cartesian(0, 45)$u, 0)
  checkEqualsNumeric(convert_polar_to_cartesian(0, 45)$v, 0)
  
  # Test for negative speed
  result_negative <- convert_polar_to_cartesian(-1, 45)
  checkEqualsNumeric(result_negative$u, -cos(45 * pi / 180), "Expected correct u for negative speed")
  checkEqualsNumeric(result_negative$v, -sin(45 * pi / 180), "Expected correct v for negative speed")
  
}

climdex.pcic.test.convert.degrees.to.cardinal <- function() {
  # Test all 16 cardinal directions and boundary cases
  degrees <- c(0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 
               202.5, 225, 247.5, 270, 292.5, 315, 337.5, 348.75)
  expected <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 
                'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  result <- convert_degrees_to_cardinal(degrees)
  
  checkEquals(result, expected)
  
  # 360-0 boundary (e.g., 360, should return 'N')
  checkEquals(convert_degrees_to_cardinal(360), 'N')
  
  # Out-of-bound degrees
  result_out_of_bound <- convert_degrees_to_cardinal(c(-10, 370))
  checkEquals(result_out_of_bound, c('N', 'N'), "Expected 'N' and 'N' for out-of-bound degrees")
}


climdex.pcic.test.convert.cardinal.to.degrees <- function() {
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

climdex.pcic.test.filter.by.direction.range <- function() {
  primary_data <- 1:10
  degrees <- seq(0, 360, length.out = 10)
  date_factors <- rep(1, 10)
  direction.range <- c(90, 270)
  
  filtered <- filter_by_direction_range(primary_data, degrees, date_factors, direction.range)
  
  expected_data <- ifelse(degrees >= 90 & degrees <= 270, primary_data, NA)
  expected_degrees <- ifelse(degrees >= 90 & degrees <= 270, degrees, NA)
  
  checkEquals(filtered$primary_data, expected_data)
  checkEquals(filtered$degrees, expected_degrees)
  
  # Test range that crosses 0 
  direction.range <- c(350, 10)  # Only include data associated with degrees 0 & 360
  filtered <- filter_by_direction_range(primary_data, degrees, date_factors, direction.range)
  
  # Expect NA for all values except for those near 0 and 360
  expected_data_cross_0 <- ifelse(degrees >= 350 | degrees <= 10, primary_data, NA)
  expected_degrees_cross_0 <- ifelse(degrees >= 350 | degrees <= 10, degrees, NA)
  
  checkEquals(filtered$primary_data, expected_data_cross_0, "Expected filtered primary data for range crossing 0 to have NA")
  checkEquals(filtered$degrees, expected_degrees_cross_0, "Expected filtered degrees for range crossing 0 to have NA")
}

climdex.pcic.test.filter.by.direction.range.full.na <- function() {
  primary_data <- 1:10
  degrees <- seq(0, 360, length.out = 10) # 0  40  80 120 160 200 240 280 320 360
  date_factors <- rep(1, 10)
  direction.range <- c(300, 310)  # Filter everything out
  
  filtered <- filter_by_direction_range(primary_data, degrees, date_factors, direction.range)
  
  expected_data <- rep(NA, 10)
  expected_degrees <- rep(NA, 10)
  
  checkEquals(filtered$primary_data, expected_data, "Expected all primary data to be NA")
  checkEquals(filtered$degrees, expected_degrees, "Expected all degrees to be NA")
}

climdex.pcic.test.compute.circular.mean <- function() {
  direction_degrees <- c(350, 10)  # Should average to 0 degrees
  date.factors <- 1
  format <- "polar"
  result <- compute_circular_mean(direction_degrees, date.factors, format)
  expected <- 0
  
  checkEqualsNumeric(result, expected, tolerance = 1e-6)
  # Single Direction
  checkEqualsNumeric(compute_circular_mean(c(90, 90), c(1, 1), "polar"), 90)
  # No directions 
  checkException(compute_circular_mean(c(), c(), "polar"), 
                 msg = "direction_degrees cannot be empty or NULL.")
}

climdex.pcic.test.compute.circular.mean.with.na <- function() {
  direction_degrees <- c(350, 10, NA)  # Should average to 0 degrees
  date.factors <- 1
  format <- "polar"
  result <- compute_circular_mean(direction_degrees, date.factors, format)
  expected <- 0
  checkEqualsNumeric(result, expected, tolerance = 1e-6)
}

climdex.pcic.test.compute.circular.sd <- function() {
  direction_degrees <- c(350, 10)
  date.factors <- 1
  result <- compute_circular_sd(direction_degrees, date.factors)
  expected <- circular::sd.circular(circular::circular(c(350, 10), units = "degrees", modulo = "2pi")) * 180 / pi
  checkEqualsNumeric(result, expected, tolerance = 1e-4)
  
  # Check for empty inputs (expect an error)
  checkException(compute_circular_sd(c(), c()), 
                 msg = "direction_degrees cannot be empty or NULL.")
}

climdex.pcic.test.compute.circular.sd.with.na <- function() {
  direction_degrees <- c(350, 10, NA)
  date.factors <- 1
  result <- compute_circular_sd(direction_degrees, date.factors)
  expected_sd <- circular::sd.circular(circular::circular(c(350, 10), units = "degrees", modulo = "2pi")) * 180 / pi
  checkEqualsNumeric(result, expected_sd, tolerance = 1e-4, "Expected circular SD ignoring NA values")
}
