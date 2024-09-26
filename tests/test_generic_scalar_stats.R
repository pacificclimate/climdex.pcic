library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"

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
