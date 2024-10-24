library(climdex.pcic)
library(RUnit)
library(circular)
# Setup
wd <- "./"

# # compute.stat.vector Tests
# Helper function to update primary and secondary values in a vector object. Dates need to be readded to avoid data-date length mismatch.
update_vector_obj <- function(base_vector_obj, primary, secondary, dates) {
  return(climdexGenericVector.raw(
    primary = primary,
    secondary = secondary,
    dates = dates,
    max.missing.days = base_vector_obj@max.missing.days,
    format = base_vector_obj@format,
    northern.hemisphere = base_vector_obj@northern.hemisphere,
    calendar = attr(base_vector_obj@dates, "cal")
  ))
}


climdex.pcic.test.compute.stat.vector.circular.mean <- function() {
  # Common setup for the vector object
  speed <- c(5, 5)  
  direction <- c(350, 10)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 2)
  
  # Create a base vector object
  base_vector_obj <- climdexGenericVector.raw(
    primary = speed,        
    secondary = direction,  
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),  
    format = "polar",       
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Test case 1: Compute circular mean for the default setup
  result <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_mean",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  expected_direction <- 0  # Circular mean of 350 and 10 degrees is 0 degrees
  checkEqualsNumeric(as.numeric(result$direction[1]), expected_direction, tolerance = 1e-6)
  
  # Test case 2: Same direction (90, 90), update the secondary component
  base_vector_obj <- update_vector_obj(base_vector_obj, primary = c(5, 5), secondary = c(90, 90), dates = dates)
  result_same_direction <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_mean",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  checkEqualsNumeric(as.numeric(result_same_direction$direction[1]), 90, tolerance = 1e-6)
  
  # Test case 3: Opposite directions (0, 180), update the secondary component
  base_vector_obj <- update_vector_obj(base_vector_obj, primary = c(5, 5), secondary = c(0, 180), dates = dates)
  result_opposite_direction <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_mean",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  checkTrue(is.na(result_opposite_direction$direction[1]), "Circular mean of opposite directions should be NA")
  
  # Test case 4: One value is NA (5, NA), update the primary and secondary components
  base_vector_obj <- update_vector_obj(base_vector_obj, primary = c(5, NA), secondary = c(350, 10), dates = dates)
  result_with_na <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_mean",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  # Should omit NA and just take mean of 350
  checkEqualsNumeric(as.numeric(result_with_na$direction[1]), 350, tolerance = 1e-6, 
                     msg = paste("Result with NA differs from expected:", result_with_na$direction[1]))
}





climdex.pcic.test.compute.stat.vector.circular.sd <- function() {
  
  speed <- c(5, 5)
  direction <- c(350, 10)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 2)
  
  base_vector_obj <- climdexGenericVector.raw(
    primary = speed,        
    secondary = direction,  
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),  
    format = "polar",       
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Test case 1: Compute circular SD for the default setup
  result <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_sd",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  
  expected_sd <-circular::sd.circular(circular::circular(c(350, 10), units = "degrees", modulo = "2pi")) * 180 / pi
  checkEqualsNumeric(as.numeric(result$circular_sd[1]), expected_sd, tolerance = 1e-4)
  
  # Test case 2: Same direction (90, 90), update the secondary component
  base_vector_obj <- update_vector_obj(base_vector_obj, primary = c(5, 5), secondary = c(90, 90), dates = dates)
  result_same_direction <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_sd",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  expected_sd_same <- 0  # SD should be 0 as both directions are the same
  checkEqualsNumeric(as.numeric(result_same_direction$circular_sd[1]), expected_sd_same, tolerance = 1e-4)
  
  # Test case 3: Opposite directions (0, 180), update the secondary component
  base_vector_obj <- update_vector_obj(base_vector_obj, primary = c(5, 5), secondary = c(0, 180), dates = dates)
  result_opposite_direction <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_sd",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  expected_sd_opp <- circular::sd.circular(circular::circular(c(0, 180), units = "degrees", modulo = "2pi")) * 180 / pi
  checkEqualsNumeric(as.numeric(result_opposite_direction$circular_sd[1]), expected_sd_opp, tolerance = 1e-6, 
                     msg = paste("Result differs from expected:", result_opposite_direction$circular_sd[1]))
  
  # Test case 4: One value is NA (5, NA), update the primary and secondary components
  base_vector_obj <- update_vector_obj(base_vector_obj, primary = c(5, NA), secondary = c(350, 10), dates = dates)
  result_with_na <- compute.stat.vector(
    base_vector_obj,
    stat = "circular_sd",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  # Should omit NA and just compute SD for the remaining value, which should be 0
  expected_sd_na <- 0
  checkEqualsNumeric(as.numeric(result_with_na$circular_sd[1]), expected_sd_na, tolerance = 1e-4)
}


climdex.pcic.test.compute.stat.vector.mean.magnitude <- function() {

  speed <- c(5, 10, 15)  
  direction <- c(0, 90, 180)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 3)

  vector_obj <- climdexGenericVector.raw(
    primary = speed,  
    secondary = direction,  
    dates = dates,  #
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),  
    format = "polar",
    northern.hemisphere = TRUE, 
    calendar = "gregorian"  
  )
  
  # Monthly mean
  result <- compute.stat.vector(
    vector_obj,
    stat = "mean",  
    freq = "monthly",  
    format = "polar",  
    include.exact.dates = FALSE
  )
  
  expected_mean_magnitude <- mean(speed)
  
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_mean_magnitude, 
                     msg = paste("Result magnitude differs from expected:", result$magnitude))
}

climdex.pcic.test.compute.stat.vector.filtered.direction.range <- function() {
  
  speed <- c(5, 10, 15, 20)
  direction <- c(45, 90, 135, 180)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 4)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Monthly max, filter by 90 to 135 degrees
  result <- compute.stat.vector(
    vector_obj,
    stat = "max",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE,
    direction.range = c(90, 135)
  )
  
  filtered_speeds <- speed[direction >= 90 & direction <= 135]
  expected_max_magnitude <- 15
  
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_max_magnitude)
}

climdex.pcic.test.compute.stat.vector.missing.values <- function() {
  speed <- c(5, 10, NA, 20) # One NA speed
  direction <- c(45, 90, 135, NA)  # One NA direction
  
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 4)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  # Monthly mean 
  result <- compute.stat.vector(
    vector_obj,
    stat = "mean",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE
  )
  valid_entries <- !is.na(speed) & !is.na(direction)
  expected_mean_magnitude <- mean(speed[valid_entries], na.rm = TRUE)  # Mean of 5 and 10 (only valid entries)
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_mean_magnitude, 
                     msg = paste("Result magnitude differs from expected. Result:", result$magnitude))
}

climdex.pcic.test.compute.stat.vector.overlapping.years <- function() {
  # Full winter season (Dec-Feb)
  set.seed(123)
  speed <- runif(90, min = 5, max = 20)  
  direction <- runif(90, min = 0, max = 360)  
  
  dates <- seq(as.PCICt("2019-12-01", cal = "gregorian"), by = "day", length.out = 90)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 3, seasonal = 6),
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Seasonal mean (expect it to compute correctly across Dec-Feb)
  result <- compute.stat.vector(
    vector_obj,
    stat = "mean",
    freq = "seasonal",
    format = "polar",
    include.exact.dates = FALSE
  )
  expected_mean_speed <- mean(speed)

  checkEqualsNumeric(as.numeric(result$magnitude["Winter 2019"]), expected_mean_speed, 
                     msg = paste("Computed seasonal mean:", result$magnitude[1], 
                                 "Expected seasonal mean:", expected_mean_speed))
}

climdex.pcic.test.compute.stat.vector.filtered.direction.crossing.360 <- function() {
  speed <- c(5, 10, 15, 20)
  direction <- c(350, 10, 135, 180)  # Directions crossing the 360-degree boundary
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 4)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Monthly max, filter by 350 to 10 degrees
  result <- compute.stat.vector(
    vector_obj,
    stat = "max",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE,
    direction.range = c(350, 10)
  )
  filtered_speeds <- speed[direction >= 350 | direction <= 10]
  expected_max_magnitude <- max(filtered_speeds)
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_max_magnitude)
}

climdex.pcic.test.compute.stat.vector.no.data.in.direction.range <- function() {
  speed <- c(5, 10, 15, 20)
  direction <- c(45, 90, 135, 180)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 4)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Monthly max, filter by 270 to 360 degrees
  result <- compute.stat.vector(
    vector_obj,
    stat = "max",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE,
    direction.range = c(270, 360)
  )

  result$magnitude[1]
  checkTrue(is.na(result$magnitude[1]), "Expected NA result when no data is in the specified direction range.")
}

climdex.pcic.test.compute.stat.vector.cartesian.format <- function() {
  u <- c(5, 10, 15)
  v <- c(5, 0, -5)  # Cartesian components
  
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 3)
  
  vector_obj <- climdexGenericVector.raw(
    primary = u,
    secondary = v,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "cartesian",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Monthly mean
  result <- compute.stat.vector(
    vector_obj,
    stat = "mean",
    freq = "monthly",
    format = "cartesian", 
    include.exact.dates = FALSE
  )
  
  # Convert u and v to magnitude and direction
  polar_data <- convert_cartesian_to_polar(u, v)
  expected_mean_magnitude <- mean(polar_data$speed)
  
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_mean_magnitude, tolerance = 1e-6,
                     msg = paste("Computed mean magnitude:", result$magnitude[1], 
                                 "Expected mean magnitude:", expected_mean_magnitude))
  # Test 2: Monthly circular mean
  result_circular_mean <- compute.stat.vector(
    vector_obj,
    stat = "circular_mean",
    freq = "monthly",
    format = "cartesian",
    include.exact.dates = FALSE
  )
  
  expected_circular_mean <- compute_circular_mean(polar_data$direction, rep(1, length(polar_data$direction)), "polar")
  
  checkEqualsNumeric(as.numeric(result_circular_mean$direction[1]), expected_circular_mean, tolerance = 1e-6,
                     msg = paste("Computed circular mean direction:", result_circular_mean$direction[1],
                                 "Expected circular mean direction:", expected_circular_mean))
}

climdex.pcic.test.compute.stat.vector.cardinal.format <- function() {
  speed <- c(5, 10, 15)
  direction <- c('N', 'E', 'S')  # Cardinal directions
  
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 3)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,        
    secondary = direction,  
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "cardinal",       
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  result <- compute.stat.vector(
    vector_obj,
    stat = "mean",
    freq = "monthly",
    format = "cardinal",
    include.exact.dates = FALSE
  )
  
  expected_mean_speed <- mean(speed)
  
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_mean_speed, tolerance = 1e-6)
}

climdex.pcic.test.compute.stat.vector.inverted.direction.range <- function() {
  speed <- c(5, 10, 15, 20)
  direction <- c(45, 90, 135, 180)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 4)
  
  vector_obj <- climdexGenericVector.raw(
    primary = speed,
    secondary = direction,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    format = "polar",
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  
  # Monthly mean, filter from 270 to 90 degrees
  result <- compute.stat.vector(
    vector_obj,
    stat = "mean",
    freq = "monthly",
    format = "polar",
    include.exact.dates = FALSE,
    direction.range = c(270, 90)
  )
  
  filtered_speeds <- speed[direction <= 90 | direction >= 270]  # Filtering with inverted bounds
  expected_mean_magnitude <- mean(filtered_speeds)
  
  checkEqualsNumeric(as.numeric(result$magnitude[1]), expected_mean_magnitude)
}
