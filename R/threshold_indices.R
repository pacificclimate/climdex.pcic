#' Frost Days
#' 
#' This function computes the climdex index FD.
#' 
#' This function takes a climdexInput object as input and computes the FD (frost
#' days) climdex index: that is, the annual count of days where daily minimum
#' temperature drops below 0 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of frost days for each year.
#' @template generic_seealso_references
#' 
#' @templateVar cdxvar fd
#' @templateVar cdxdescription an annual timeseries of the number of frost days.
#' @template get_generic_example
#' 
#' @export
climdex.fd <- function(ci) { stopifnot(!is.null(ci@data$tmin)); return(number.days.op.threshold(ci@data$tmin, ci@date.factors$annual, 0, "<") * ci@namasks$annual$tmin) }

#' Summer Days
#' 
#' This function computes the climdex index SU.
#' 
#' This function takes a climdexInput object as input and computes the SU (summer
#' days) climdex index: that is, the annual count of days where daily maximum
#' temperature exceeds 25 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of summer days for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar su
#' @templateVar cdxdescription an annual timeseries of the number of summer days.
#' @template get_generic_example
#' 
#' @export
climdex.su <- function(ci) { stopifnot(!is.null(ci@data$tmax)); return(number.days.op.threshold(ci@data$tmax, ci@date.factors$annual, 25, ">") * ci@namasks$annual$tmax) }

#' Icing Days
#' 
#' This function computes the climdex index ID.
#' 
#' This function takes a climdexInput object as input and computes the ID (icing
#' days) climdex index: that is, the annual count of days where daily maximum
#' temperature is below 0 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of icing days for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar id
#' @templateVar cdxdescription an annual timeseries of the number of icing days.
#' @template get_generic_example
#' 
#' @export
climdex.id <- function(ci) { stopifnot(!is.null(ci@data$tmax)); return(number.days.op.threshold(ci@data$tmax, ci@date.factors$annual, 0, "<") * ci@namasks$annual$tmax) }

#' Tropical Nights
#' 
#' This function computes the climdex index TR.
#' 
#' This function takes a climdexInput object as input and computes the TR
#' (tropical nights) climdex index: that is, the annual count of days where
#' daily minimum temperature stays above 20 degrees Celsius.
#' 
#' @param ci Object of type climdexInput.
#' @return A vector containing the number of tropical nights for each year.
#' @template generic_seealso_references
#' @templateVar cdxvar tr
#' @templateVar cdxdescription an annual timeseries of the number of tropical nights.
#' @template get_generic_example
#' 
#' @export
climdex.tr <- function(ci) { stopifnot(!is.null(ci@data$tmin)); return(number.days.op.threshold(ci@data$tmin, ci@date.factors$annual, 20, ">") * ci@namasks$annual$tmin) }

#' @title Growing Season Length
#' 
#' @description
#' This function computes the growing season length (GSL) given the input.
#' 
#' @details
#' This function takes a climdexInput object as input and computes the growing
#' season length based on this data.
#' 
#' Growing season length as defined by the climdex indices is the number of
#' days between the start of the first spell of warm days in the first half of
#' the year, and the start of the first spell of cold days in the second half
#' of the year. Spells of warm days are defined as six or more days with mean
#' temperature above 5 degrees Celsius; spells of cold days are defined as six
#' or more days with a mean temperature below 5 degrees Celsius.
#' 
#' The three alternate modes provided ('GSL_first', 'GSL_max', and 'GSL_sum')
#' are for testing purposes only. They differ considerably from the first
#' ('GSL') mode. All of them use a list of growing seasons -- here defined as
#' six or more consecutive days with a mean temperature greater than or equal
#' to 5 degrees Celsius, followed by either the end of the year or six or more
#' consecutive days with a mean temperature less than 5 degrees Celsius.
#' 'GSL_first' returns the first growing season found; 'GSL_max' returns the
#' longest growing season found; and 'GSL_sum' returns the total length of all
#' growing seasons found.
#' 
#' @param ci Object of type climdexInput.
#' @param gsl.mode Growing season length method to use.
#' @param include.exact.dates Logical, if TRUE, return a data frame with the season lengths and the start and end dates of each season; if FALSE, return only the season lengths.
#' @return A vector containing the number of days in the growing season for
#' each year.
#' @note Note that fclimdex results may differ from results using the first
#' ('GSL') mode due to bugs in fclimdex. Please ensure you are using the latest
#' version of fclimdex, as there have been numerous bug fixes and the results
#' should, at this point, match.
#' 
#' Please do not use the 'GSL_first', 'GSL_max', or 'GSL_sum' modes for
#' anything other than testing purposes at this time, nor should you rely on
#' this parameter being present in future versions of climdex.pcic.
#' @seealso \code{\link{growing.season.length}},
#' \code{\link{climdexInput.csv}}.
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @keywords ts climate
#' @templateVar cdxvar gsl
#' @templateVar cdxdescription an annual timeseries of the growing season length in days.
#' @template get_generic_example
#' 
#' @export

climdex.gsl <- function(ci, gsl.mode = c("GSL", "GSL_first", "GSL_max", "GSL_sum"), include.exact.dates = FALSE) {
  stopifnot(!is.null(ci@data$tavg))
  ## Gotta shift dates so that July 1 is considered Jan 1 of same year in southern hemisphere
  if (ci@northern.hemisphere) {
    gs <- growing.season.length(ci@data$tavg, ci@date.factors$annual, ci@dates, ci@northern.hemisphere, gsl.mode = match.arg(gsl.mode), include.exact.dates = include.exact.dates, cal = attr(ci@dates, "cal"))
    namask.gsl <- ci@namasks$annual$tavg
    
  } else {
    dates.POSIXlt <- as.POSIXlt(ci@dates)
    years <- dates.POSIXlt$year + 1900
    months <- dates.POSIXlt$mon + 1
    
    valid.years <- range(years)
    years.gsl <- years - floor((12 - months) / 6)
    
    inset <- years.gsl >= valid.years[1]
    gsl.factor <- factor(years.gsl[inset])
    gsl.factor.monthly <- factor(paste(years.gsl[inset], months[inset], sep = "-"))
    gsl.yearmonth.factor <- unlist(strsplit(levels(gsl.factor.monthly), "-"))[(0:(nlevels(gsl.factor.monthly) - 1)) * 2 + 1]
    gsl.temp.data <- ci@data$tavg[inset]
    namask.gsl.monthly <- get.na.mask(gsl.temp.data, gsl.factor.monthly, ci@max.missing.days["annual"])
    namask.gsl <- get.na.mask(gsl.temp.data, gsl.factor, ci@max.missing.days["annual"]) * as.numeric(tapply(namask.gsl.monthly, gsl.yearmonth.factor, prod))
    dim(namask.gsl) <- dimnames(namask.gsl) <- NULL
    namask.gsl[length(namask.gsl)] <- NA
    
    gs <- growing.season.length(gsl.temp.data, gsl.factor, ci@dates[inset], ci@northern.hemisphere, gsl.mode = match.arg(gsl.mode), include.exact.dates = include.exact.dates, cal = attr(ci@dates, "cal"))
  }
  if (include.exact.dates) {
    gs$sl <- gs$sl * namask.gsl
    gs$start[is.na(gs$sl)] <- NA
    gs$end[is.na(gs$sl)] <- NA
    return(gs)
  }
  
  return(gs * namask.gsl)
}

#' @title Warm Spell Duration Index
#'
#' @description This function computes the climdex index WSDI.
#' 
#' @details
#' This function takes a climdexInput object as input and computes the climdex
#' index WSDI (Warm Spell Duration Index).
#' 
#' The warm spell duration index is defined as the number of days each year
#' which are part of a "warm spell". A "warm spell" is defined as a sequence of
#' 6 or more days in which the daily maximum temperature exceeds the 90th
#' percentile of daily maximum temperature for a 5-day running window
#' surrounding this day during the baseline period.
#' 
#' The \code{spells.can.span.years} option specifies whether spells can cross
#' year boundaries -- i.e., span years. The default for this is the same as
#' fclimdex.
#' 
#' @template wcsdi_common
#' @templateVar cdxvar wsdi
#' @templateVar cdxdescription an annual timeseries of the warm spell duration index.
#' @template get_generic_example
#' 
#' @export
climdex.wsdi <- function(ci, spells.can.span.years=FALSE) { stopifnot(!is.null(ci@data$tmax) && !is.null(ci@quantiles$tmax)); return(threshold.exceedance.duration.index(ci@data$tmax, ci@date.factors$annual, ci@jdays, ci@quantiles$tmax$outbase$q90, ">", spells.can.span.years=spells.can.span.years, max.missing.days=ci@max.missing.days['annual']) * ci@namasks$annual$tmax) }

#' @title Cold Spell Duration Index
#' 
#' @description This function computes the climdex index CSDI.
#' 
#' @details
#' This function takes a climdexInput object as input and computes the climdex
#' index CSDI (Cold Spell Duration Index).
#'
#' The cold spell duration index is defined as the number of days
#' each year which are part of a "cold spell". A "cold spell" is defined as a
#' sequence of 6 or more days in which the daily minimum temperature is below
#' the 10th percentile of daily minimum temperature for a 5-day running window
#' surrounding this day during the baseline period.
#' 
#' The \code{spells.can.span.years} option specifies whether spells can cross
#' year boundaries -- i.e., span years. The default for this is the same as
#' fclimdex.
#' 
#' @template wcsdi_common
#' @templateVar cdxvar csdi
#' @templateVar cdxdescription an annual timeseries of the cold spell duration index.
#' @template get_generic_example
#' 
#' @export
climdex.csdi <- function(ci, spells.can.span.years=FALSE) { stopifnot(!is.null(ci@data$tmin) && !is.null(ci@quantiles$tmin)); return(threshold.exceedance.duration.index(ci@data$tmin, ci@date.factors$annual, ci@jdays, ci@quantiles$tmin$outbase$q10, "<", spells.can.span.years=spells.can.span.years, max.missing.days=ci@max.missing.days['annual']) * ci@namasks$annual$tmin) }


#' @title Flexible GSL function
#' 
#' @description
#' This function computes the growing season length (GSL) given the input,
#' which is allowed to vary considerably from the ETCCDI definitions.
#' 
#' @details
#' This function is the function used to implement \code{\link{climdex.gsl}}.
#' It's designed to be flexible to allow for experimentation and testing of new
#' thresholds and methods.
#' 
#' If you need to use this code for experimentation in the southern hemisphere,
#' you'll need to rip off the climdex.gsl code to rotate the year around so
#' that July 1st is treated as January 1st.
#' 
#' See \code{\link{climdex.gsl}} for more information on what \code{gsl.mode}
#' does.
#' 
#' @param daily.mean.temp Timeseries of daily mean temperature (in degrees C),
#' padded out to end on a year boundary (ie: starts on January 1st of some
#' year, ends on December 31st).
#' @param date.factor Factor of the same length as daily.mean.temp that divides
#' the timeseries up into years of data.
#' @param dates The corresponding series of dates.
#' @param northern.hemisphere Whether the data is from the northern hemisphere.
#' @param min.length The minimum number of days above or below the threshold
#' temperature that defines the start or end of a growing season.
#' @param t.thresh The temperature threshold for being considered part of a
#' growing season (in degrees C).
#' @param gsl.mode The growing season length mode (ETCCDI mode is "GSL").
#' @param include.exact.dates Logical, if TRUE, return a data frame with the season lengths and the start and end dates of each season; if FALSE, return only the season lengths.
#' @param cal Calendar object specifying the calendar system.
#' @return A vector containing the number of days in the growing season for
#' each year.
#' @seealso \code{\link{climdex.gsl}}, \code{\link{climdexInput.csv}}.
#' @keywords ts climate
#' @examples
#' library(PCICt)
#' 
#' ## Create a climdexInput object from some data already loaded in and
#' ## ready to go.
#' 
#' ## Parse the dates into PCICt.
#' tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' prec.dates <- as.PCICt(do.call(paste, ec.1018935.prec[,c("year",
#' "jday")]), format="%Y %j", cal="gregorian")
#' 
#' ## Load the data in.
#' ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP,
#' ec.1018935.tmin$MIN_TEMP, ec.1018935.prec$ONE_DAY_PRECIPITATION,
#' tmax.dates, tmin.dates, prec.dates, base.range=c(1971, 2000))
#' 
#' ## Create an annual timeseries of the growing season length in days.
#' gsl <- growing.season.length(ci@@data$tavg, ci@@date.factors$annual, ci@@dates,
#'                              ci@@northern.hemisphere, gsl.mode="GSL") * 
#'        ci@@namasks$annual$tavg
#' 
#' ## Print these out for testing purposes.
#' gsl
#' 
#' @export
growing.season.length <- function(daily.mean.temp, date.factor, dates, northern.hemisphere,
                                  min.length = 6, t.thresh = 5, gsl.mode = c("GSL", "GSL_first", "GSL_max", "GSL_sum"), include.exact.dates = FALSE, cal) {
  gsl.mode <- match.arg(gsl.mode)
  month.series <- get.months(dates)
  transition.month <- if (northern.hemisphere) 7 else 1
  if (gsl.mode == "GSL") {
    growing.season <- (tapply.fast(1:length(daily.mean.temp), date.factor, function(idx) {
      temp.data <- daily.mean.temp[idx]
      ts.mid <- head(which(month.series[idx] == transition.month), n = 1)
      if (!length(ts.mid)) {
        return(NA)
      }
      
      ts.len <- length(temp.data)
      gs.begin <- which(select.blocks.gt.length(temp.data[1:(ts.mid - 1)] > t.thresh, min.length - 1))
      
      ## Growing season actually ends the day -before- the sequence of sketchy days
      gs.end <- which(select.blocks.gt.length(temp.data[ts.mid:ts.len] < t.thresh, min.length - 1)) - 1
      
      ## If no growing season start, 0 length; if no end, ends at end of year; otherwise, end - start + 1
      # season.length <- ifelse(length(gs.begin) == 0, 0, ifelse(length(gs.end) == 0, ts.len - gs.begin[1] + 1, gs.end[1] - gs.begin[1] + ts.mid))
      
      if (is.na(gs.begin[1])) {
        start <- NA_character_
        end <- NA_character_
        season.length <- 0  # Set season length to 0
      } 
      else if (is.na(gs.end[1])) {
        start <- gs.begin[1]
        end <- ts.len - 1 # Last DOY
        season.length <- ts.len - start + 1
      } 
      else {
        start <- gs.begin[1]
        end <- gs.end[1] + ts.mid -1
        season.length <- end - start +1 
      }
      
      if (include.exact.dates) {
        if(northern.hemisphere){
          origin <- paste(year = unique(date.factor[idx]), "01", "01", sep = "-")
        }
        else {origin <- paste(year = unique(date.factor[idx]), "07", "01", sep = "-")}
        
        origin.pcict <- as.PCICt(origin, cal)
        seconds.per.day <- 86400
        
        if (!is.na(start)){
          start.pcict <- origin.pcict + (start - 1) * seconds.per.day
          start <- format(start.pcict, "%Y-%m-%d")
          end.pcict <- origin.pcict + (end * seconds.per.day)
          end <- format(end.pcict, "%Y-%m-%d")
        }
        
        result <- c(start, sl <- season.length, end)
      } 
      else { result <- season.length}
      return(result)
    }))
    
    if (include.exact.dates) {
      start <- growing.season[seq(1, length(growing.season), by = 3)]
      sl <- growing.season[seq(2, length(growing.season), by = 3)]
      sl <- as.numeric(sl)
      end <- growing.season[seq(3, length(growing.season), by = 3)]
      year <- names(growing.season[1:length(start)]) 
      if(!northern.hemisphere){
        sl <- c(sl,NA)
        end <- c(end,NA)
      }
      df <- data.frame(start,sl, end)
      rownames(df) <- year
      return(df)
    }
    return(growing.season)
  } 
  else {
    in.gsl <- !select.blocks.gt.length(!select.blocks.gt.length(daily.mean.temp >= t.thresh, min.length - 1), min.length - 1)
    warning("GSL_first, GSL_max, and GSL_sum are experimental alternative growing season length definitions. Use at your own risk.")
    
    innerfunc <- switch(gsl.mode,
                        GSL_first = function(bl) {
                          ifelse(any(bl > 0), (bl[bl > 0])[1], 0)
                        },
                        GSL_max = max,
                        GSL_sum = sum
    )
    return(tapply.fast(in.gsl, date.factor, function(ts) {
      block.lengths <- get.series.lengths.at.ends(ts)
      return(innerfunc(block.lengths))
    }))
  }
}
