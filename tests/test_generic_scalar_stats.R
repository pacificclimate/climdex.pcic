library(climdex.pcic)
library(RUnit)

# Setup
wd <- "./"

# # compute.stat.scalar Tests

climdex.pcic.test.compute.stat.scalar.mean <- function() {
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
  # Monthly mean
  result <- compute.stat.scalar(scalar_obj, stat = "mean", freq = "monthly", include.exact.dates = FALSE)
  expected_mean <- mean(data)
  
  checkEqualsNumeric(as.numeric(result[1]), expected_mean)
}


climdex.pcic.test.compute.stat.scalar.max <- function() {
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
  
  # Monthly max
  result <- compute.stat.scalar(scalar_obj, stat = "max", freq = "monthly", include.exact.dates = FALSE)
  expected_max <- max(data)
  
  checkEqualsNumeric(as.numeric(result[1]), expected_max)
}

climdex.pcic.test.scalar.exact.dates <- function() {
  set.seed(123)
  
  data <- runif(365, 0, 20)
  dates <- seq(as.PCICt("2020-01-01", cal = "gregorian"), by = "day", length.out = 365)
  
  scalar_obj <- climdexGenericScalar.raw(
    data = data,
    dates = dates,
    max.missing.days = c(annual = 15, monthly = 31, seasonal = 6),
    northern.hemisphere = TRUE,
    calendar = "gregorian"
  )
  # Monthly max
  result <- compute.stat.scalar(scalar_obj, stat = "max", freq = "monthly", include.exact.dates = TRUE)
  
  expected_max <- max(data)
  max_index <- which.max(data)
  expected_max_date <- dates[max_index]
  
  computed_exact_date <- as.PCICt(result$ymd[1], cal = "gregorian")
  
  checkEqualsNumeric(as.numeric(result$val[1]), expected_max, 
                     msg = paste("Computed max statistic:", result$val[1], 
                                 "Expected max statistic:", expected_max))
  
  checkEquals(computed_exact_date, expected_max_date, 
              msg = paste("Computed exact date:", computed_exact_date, 
                          "Expected exact date:", expected_max_date))
}

