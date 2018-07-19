library(climdex.pcic)
library(RUnit)

climdex.pcic.test.running.mean <- function() {
  vec <- 1:10
  error.message <- "running.mean failed with window of width "
  #handle window width 1
  means <- climdex.pcic:::running.mean(vec, 1)
  checkEquals(means, 1:10, paste(error.message, 1))
  
  #handle odd and even widths
  means <- climdex.pcic:::running.mean(vec, 2)
  checkEquals(means, c(NA, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5), paste(error.message, 2))
  means <- climdex.pcic:::running.mean(vec, 3)
  checkEquals(means, c(NA, 2, 3, 4, 5, 6, 7, 8, 9, NA), paste(error.message, 3))
  
  #handle width equal to or greater than the data
  means <- climdex.pcic:::running.mean(vec, 10)
  checkEquals(means, c(NA, NA, NA, NA, NA, 5.5, NA, NA, NA, NA), paste(error.message, "maximum"))
  means <- climdex.pcic:::running.mean(vec, 11)
  checkEquals(means, climdex.pcic:::running.mean(vec, 10), paste(error.message, "exceeding maximum"))
}