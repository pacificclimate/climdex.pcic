library(climdex.pcic)
library(RUnit)
library(circular)
# Setup
wd <- "./"

# 
# # compute.stat.vector Tests

# Helper function to update primary and secondary values in a vector object. Dates need to be readded to avoid data-date length 
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


climdex.pcic.test.compute_stat_vector_circular_mean <- function() {
  # Common setup for the vector object
  speed <- c(5, 5)  # Speed values
  direction <- c(350, 10)  # Direction values in degrees
  
  # Use PCICt and seq to create a date sequence
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 2)
  
  # Create a base vector object using ClimdexGenericVector.raw
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





climdex.pcic.test.compute_stat_vector_circular_sd <- function() {
  # Common setup for the vector object
  speed <- c(5, 5)  # Speed values
  direction <- c(350, 10)  # Direction values in degrees
  
  # Use PCICt and seq to create a date sequence
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 2)
  
  # Create a base vector object using ClimdexGenericVector.raw
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


climdex.pcic.test.compute_stat_vector_mean_magnitude <- function() {
  # Test data: speed and direction values
  speed <- c(5, 10, 15)  
  direction <- c(0, 90, 180)
  
  # Generate dates using PCICt (calendar type: gregorian)
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
  
  # Compute the mean magnitude for the 'monthly' frequency
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

climdex.pcic.test.compute_stat_vector_filtered_direction_range <- function() {
  
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
  
  # Compute the max statistic, filtered by the direction range (90 to 135 degrees)
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


