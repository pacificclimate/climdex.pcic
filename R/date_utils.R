#' Get the last month and day of the year
#'
#' Get the last month and day of the year as a character sting, separated by
#' the specified separator.
#'
#' This is a utility function necessitated by 360-day calendars. Works on PCICt objects.
#'
#' @param d An exemplar date.
#' @param sep Separator to use.
#' @return A string (like "12-30", or "12-31")
#' 
#' @examples
#' library(PCICt)
#' last.mday <- get.last.monthday.of.year(as.PCICt("2011-01-01", cal="360"))
#' 
#' @export
get.last.monthday.of.year <- function(d, sep="-") {
  if(!is.null(attr(d, "months"))) paste("12", attr(d, "months")[12], sep=sep) else paste("12", "31", sep=sep)
}

## Get julian day of year
get.jdays <- function(dates) {
  return(as.POSIXlt(dates)$yday + 1)
}

## Get year
get.years <- function(dates) {
  return(as.POSIXlt(dates)$year + 1900)
}

## Get month number
get.months <- function(dates) {
  return(as.POSIXlt(dates)$mon + 1)
}

## Juggle the list so that day 366 == day 365
get.jdays.replaced.feb29 <- function(jdays) {
  indices <- which(jdays == 366)
  if(length(indices) > 0)
    jdays[rep(indices, each=366) + -365:0] <- c(1:59, 59, 60:365)
  jdays
}

# Converts a positional index with respect to some origin into a PCICt object in the format %Y-%m-%d.
ymd.dates <- function(origin, cal, exact.day, val) {
  origin.pcict <- as.PCICt(origin, cal)
  seconds.per.day <- 86400
  exact.day.pcict <- origin.pcict + (ifelse(is.na(exact.day),1,exact.day - 1)) * seconds.per.day
  ymd <- as.PCICt(exact.day.pcict, cal = cal, format = "%Y-%m-%d")
  ymd <- format(ymd, "%Y-%m-%d")
  ymd[is.na(val)] <- NA
  return(ymd)
}


# Computes exact dates for statistics based on the specified frequency (annual, monthly, or seasonal).
exact.date <- function(stat, data, date.factor, freq, cal, mask) {
  val <- suppressWarnings(tapply.fast(data, date.factor, get(stat), na.rm = TRUE)) * mask
  exact.day <- suppressWarnings(tapply.fast(data, date.factor, get(paste("which.", stat, sep = "")))) 
  df <- data.frame(
    val = val,
    ymd = {
      origin <- sapply(1:length(unique(date.factor)), function(i) {
        switch(
          as.character(freq),
          annual = paste((unique(date.factor))[[i]], "01-01", sep = "-"),
          monthly = paste((unique(date.factor))[[i]], "01", sep = "-"),
          seasonal = {
            season.year <- strsplit(as.character(unique(date.factor)[[i]]), " ")
            year <- as.numeric(season.year[[1]][2])
            season <- season.year[[1]][1]
            season.months <- list(Winter = "12", Spring = "03", Summer = "06", Fall = "09")
            paste(year, season.months[[season]], "01", sep = "-")
          }
        )
      })
      ymd.dates(origin, cal, exact.day, val)
    }
  )
  return(df)
}

# Classify meteorological seasons with associated years. This includes 
# handling for the winter season, where data in the months of January and 
# February are assigned to the winter season of the previous year.
# 
# Seasons are defined as the meteorological seasons:
#  - 'Winter': December, January, February
#  - 'Spring': March, April, May
#  - 'Summer': June, July, August
#  - 'Fall': September, October, November
classify_meteorological_season_with_year <- function(dates) {
  month <- as.integer(format(dates, "%m"))
  year <- as.integer(format(dates, "%Y"))
  year_for_season <- ifelse(month %in% c(1, 2), year - 1, year)
  
  season <- ifelse(month %in% c(12, 1, 2), "Winter",
                   ifelse(month %in% 3:5, "Spring",
                          ifelse(month %in% 6:8, "Summer",
                                 ifelse(month %in% 9:11, "Fall", NA))))
  
  season_with_year <- paste(season, year_for_season)
  return(season_with_year)
}