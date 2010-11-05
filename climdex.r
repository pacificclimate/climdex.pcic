library(caTools)

## Temperature units: degrees C
## Precipitation units: mm per unit time

## FD, ID
number.days.below.threshold <- function(temp, date.factor, threshold) {
  return(tapply(temp < threshold, date.factor, sum))
}

## SU, TR
number.days.over.threshold <- function(temp, date.factor, threshold) {
  return(tapply(temp > threshold, date.factor, sum))
}

## GSL
## Meaningless if not annual
growing.season.length <- function(daily.mean.temp, date.factor, min.length=6) {
  return(tapply(daily.mean.temp, date.factor, function(ts) {
    ts.len<- length(ts)
    ts.mid <- floor(ts.len / 2)
    gs.begin <- which(select.blocks.gt.length(ts > 5, min.length - 1))
    gs.end <- which(select.blocks.gt.length(ts[ts.mid:ts.len] < 5, min.length - 1))
    #browser()
    if(length(gs.begin) == 0) {
      return(0)
    } else if(length(gs.end) == 0) {
      return(ts.len - gs.begin[1] + 1)
    } else {
      return(gs.end[1] - gs.begin[1] + 1 + ts.mid)
    }
  } ))
}

## TNx, TXx
max.daily.temp <- function(daily.temp, date.factor) {
  return(tapply(daily.temp, date.factor, max))
}

## TNn, TXn
min.daily.temp <- function(daily.temp, date.factor) {
  return(tapply(daily.temp, date.factor, min))
}

## TN10p, TX10p
## Requires use of bootstrap procedure to generate 1961-1990 pctile; see Zhang et al, 2004
percent.days.lt.threshold <- function(temp, date.factor, threshold) {
  return(tapply(temp < threshold, date.factor, function(x) { return(sum(x) / length(x) * 100) } ))
}

## TN90p, TX90p
## Requires use of bootstrap procedure to generate 1961-1990 pctile; see Zhang et al, 2004
percent.days.gt.threshold <- function(temp, date.factor, threshold) {
  return(tapply(temp > threshold, date.factor, function(x) { return(sum(x) / length(x) * 100) } ))
}

## WSDI
## Thresholds appear to be for each block of 5 days of a year...
warm.spell.duration.index <- function(daily.max.temp, dates, date.factor, warm.thresholds, min.length=6) {
  jday <- as.POSIXlt(dates)$yday + 1
  warm.periods <- select.blocks.gt.length(daily.max.temp > warm.thresholds[jday], min.length - 1)
  return(tapply(warm.periods, date.factor, sum))
}

## CSDI
## Thresholds appear to be for each block of 5 days of a year...
cold.spell.duration.index <- function(daily.min.temp, dates, date.factor, cold.thresholds) {
  jday <- as.POSIXlt(dates)$yday + 1
  cold.periods <- select.blocks.gt.length(daily.max.temp < cold.thresholds[jday], min.length - 1)
  return(tapply(cold.periods, date.factor, sum))
}

## DTR
## Max and min temps are assumed to be same length
mean.daily.temp.range <- function(daily.max.temp, daily.min.temp, date.factor) {
  return(tapply(daily.max.temp - daily.min.temp, date.factor, mean))
}

## Rx1day, Rx5day
max.nday.consec.prec <- function(daily.prec, date.factor, ndays) {
  if(ndays == 1) {
    return(tapply(daily.prec, date.factor, max))
  } else {
    ## Ends of the data will be de-emphasized (padded with zero precip data)
    prec.runsum <- runmean(c(rep(0, floor(ndays / 2)), daily.prec, rep(0, floor(ndays / 2))), k=ndays, endrule="trim") * ndays
    return(tapply(prec.runsum, date.factor, max))
  }
}

## SDII
simple.precipitation.intensity.index <- function(daily.prec, date.factor) {
  return(tapply(daily.prec, date.factor, function(prec) { idx <- prec >= 1; return(sum(prec[idx]) / sum(idx)) } ))
}

## R10mm, R20mm, Rnnmm
count.days.ge.threshold <- function(daily.prec, date.factor, threshold) {
  return(tapply(daily.prec >= threshold, date.factor, sum))
}

## CDD
max.length.dry.spell <- function(daily.prec, date.factor) {
  return(tapply(daily.prec < 1, date.factor, function(x) { return(max(sequential(x))) } ))
}

## CWD
max.length.wet.spell <- function(daily.prec, date.factor) {
  return(tapply(daily.prec >= 1, date.factor, function(x) { return(max(sequential(x))) } ))
}

## R95pTOT, R99pTOT
total.precip.above.threshold <- function(daily.prec, date.factor, threshold) {
  return(tapply(daily.prec, date.factor, function(x) { return(sum(daily.prec[daily.prec > threshold])) } ))
}

## PRCPTOT
total.prec <- function(daily.prec, date.factor) {
  return(tapply(daily.prec, date.factor, sum))
}

## Takes a list of booleans; returns a list of booleans where only blocks of TRUE longer than n are still TRUE
select.blocks.gt.length <- function(d, n) {
  if(n == 0)
    return(d)

  if(n >= length(d))
    return(rep(FALSE, length(d)))

  d2 <- Reduce(function(x, y) { return(c(rep(0, y), d[1:(length(d) - y)]) & x) }, 1:n, d)
  return(Reduce(function(x, y) { return(c(d2[(y + 1):length(d2)], rep(0, y)) | x) }, 1:n, d2))
}

## Input vector of booleans
## Returns a vector of integers representing the _length_ of each consecutive sequence of True values
sequential <- function(v) {
  if (! any(v, na.rm=T)) return(0)
  vect <- which(v)
  diff(which(c(T, diff(vect) != 1, T)))
}

## Input vector of booleans
## Returns the indicies for sequences of True which are greater than or equal to length "len"
sequential.return.indicies <- function(v, len=7) {
  i <- which(v)
  ## Represents the length-1 of each sequence of repeat meausurements
  s <- sequential(v)
  ## Which sequenses match our length critereon
  matches <- which(s >= len)

  if (length(matches) == 0) {
    return(NULL)
  }
  ## Use this as an index into our indices of repeat measurements
  match.i <- sapply(matches, function (m) {sum(s[1:m-1]) + 1})
  lengths <- s[matches]
  ## And finally return the indicies which correspond to the sequenses
  unlist(mapply(seq, i[match.i], i[match.i] + lengths - 1, SIMPLIFY=F ))
}

days.with.max.temp <- function(data, max.temp=35.0, frequency=1.0) {
  return (length(which(data > max.temp)) / frequency)
}

days.with.min.temp <- function(data, min.temp=-30.0, frequency=1.0) {
  return (length(which(data < min.temp)) / frequency)
}

long.hot.period <- function(data, max.temp=30.0, duration=7, return.only.num.events=F, return.event.indicies=F) {
  hot.days <- data > max.temp
  sequenses <- sequential(hot.days)
  sequenses <- sequenses[which(sequenses > duration)]
  if (return.only.num.events) {
    return(length(sequenses))
  }
  i <- sequential.return.indicies(hot.days, duration)
  return(c(num.events=length(sequenses), mean.magnitude=mean(sequenses, na.rm=T), indicies=i))
}

long.hot.period.indicies <- function(data, max.temp=30.0, duration=7) {
  hot.days <- data > max.temp
  sequenses <- sequential(hot.days)
  sequenses <- sequenses[which(sequenses > duration)]
  sequential.return.indicies(hot.days, duration)
}

daily.variation <- function(min.obs, max.obs, temp.range=25.0, frequency=1.0) {
  if (any(max.obs < min.obs))
    warning("Some of the min temperatures where greater than the max temperatures!")

  return (length(which(abs(max.obs - min.obs) > temp.range)) / frequency)
}

freeze.thaw <- function(min.obs, max.obs, freeze.temp=0.0, duration=85.0) {
  freeze.and.thaw <- min.obs < freeze.temp & max.obs > freeze.temp
  sequenses <- sequential(freeze.and.thaw)
  sequenses <- sequenses[which(sequenses > duration)]
  return(c(num.events=length(sequenses), mean.magnitude=mean(sequenses)))
}

deep.freeze <- function(min.obs, freeze.temp=0.0, duration=47.0) {
  frozen <- min.obs < freeze.temp
  sequenses <- sequential(frozen)
  sequenses <- sequenses[which(sequenses > duration)]
  return(c(num.events=length(sequenses), mean.magnitude=mean(sequenses)))
}

## Set duration to NA if you just want to know how many obs are greater than the threshold (i.e. no sequenses)
big.rain <- function(pcp, duration=5.0, pcp.thresh=25.0, freeze.temp=0.0) {
  rainy.days <- pcp > pcp.thresh
  if (is.na(duration)) {
    return(length(which(rainy.days)))
  }
  sequenses <- sequential(rainy.days)
  sequenses <- sequenses[which(sequenses > duration)]
  return(c(num.events=length(sequenses), mean.magnitude=mean(sequenses)))
}

blizzard <- function(pcp, temp, wind, pcp.thresh = 25.0, freeze.temp=0.0, wind.thresh=10.0, duration=24) {
  blizzard.obs <- pcp > pcp.thresh & temp < freeze.temp & wind > wind.thresh
  sequenses <- sequential(blizzard.obs)
  sequenses <- sequenses[which(sequenses > duration)]
  return(c(num.events=length(sequenses), mean.magnitude=mean(sequenses)))
}

snow <- function(pcp, temp, pcp.thresh = 10.0, freeze.temp=0.0, frequency=1.0) {
  if (length(pcp) != length(temp)) {
    warning(paste("pcp length =", length(pcp), "while temp length =", length(temp), "\nI'm going to clip them"))
    min.length <- min(length(pcp), length(temp))
    pcp <- pcp[1:min.length]
    temp <- temp[1:min.length]
  }
  snow <- pcp > pcp.thresh & temp < freeze.temp
  return(length(which(snow) / frequency))
}

pinapple.express <- function(u, v, pcp, temp, wind.speed, pcp.thresh = 25.0, temp.thresh=0.0, wind.speed.thresh, duration=24) {
  lengths <- sapply(list(u, v, pcp, temp, wind.speed), length)

  if (! all(diff(lengths) == 0)) { # lengths should be the same
    warning(paste(c("lengths of the vectors differ", lengths, "\nI'm going to clip them to the shortest length"), sep=" "))
    n <- min(lengths)
    u    <- u[1:n]
    v    <- v[1:n]
    pcp  <- pcp[1:n]
    temp <- temp[1:n]
    wind.speed <- wind.speed[1:n]
  }
  
  ## The pineapple express is defined as being from the southwest.  u and v must both be in the positive quadrant
  pina.is.coming <- u > 0 & v > 0 & temp > temp.thresh & pcp > pcp.thresh & wind.speed > wind.speed.thresh

  sequences <- sequential(pina.is.coming)
  sequences <- sequences[which(sequences > duration)]
  return(c(num.events=length(sequences), mean.magnitude=mean(sequences)))
}

rain.on.frozen.ground <- function(pcp, ts, snd, pcp.thresh=29.8, freeze.temp=0, return.indicies=F) {
  d.snd <- c(0, diff(snd))
  not.snowing <- d.snd <= 0
  hits <- pcp > pcp.thresh & ts < freeze.temp & not.snowing
  if (return.indicies) {
    return(which(hits))
  }
  else {
    return(length(which(hits)))
  }
}

rapid.snow.melt <- function(snm, snm.thresh=10.0, return.indicies=F) {
  hits <- snm > snm.thresh
  if (return.indicies) {
    return(which(hits))
  }
  else {
    return(length(which(hits)))
  }
}

total.annual.rainfall <- function(pcp, years.factor) {
  tapply(pcp, years.factor, sum)
}

high.wind <- function(w, wind.thresh=55) {
  hits <- w > wind.thresh
  return(length(which(hits)))
}

high.wind.new.dir <- function(w, dirs, wind.thresh=65, old.dir="N") {
  hits <- (w > wind.thresh) & (dirs != old.dir)
  return(length(which(hits)))
}

## If tas (air temperature) is provided, check to make sure that it is above freezing... i.e. the precip is rain and not snow/sleet/hail
rain.on.snow <- function(pcp, snd, tas=NULL, pcp.thresh=29.8, snd.thresh=10, freeze.temp=0, return.indicies=F) {
  if (is.null(tas))
    hits <- pcp > pcp.thresh & snd > snd.thresh
  else
    hits <- pcp > pcp.thresh & snd > snd.thresh & tas > freeze.temp

  if (return.indicies) {
    return(which(hits))
  }
  else {
    return(length(which(hits)))
  }
}

## FIXME: percentiles is not really a necessary argument anymore
## but it does have the same dimensionality which we want for the return value (minus the variable dimension)
## Maybe just pass in the dimnames?
apply.pina <- function(percentiles, data.list) {

  i <- which(names(dimnames(percentiles)) == "var")
  d <- dim(percentiles)[-i]
  results <- rep(NA, prod(d))
  dim(results) <- d
  dimnames(results) <- dimnames(percentiles)[-i]

  for (p in dimnames(percentiles)$percentiles) {
    for (rcm in dimnames(percentiles)$rcm) {
      for (season in dimnames(percentiles)$season) {
        for (gcm in dimnames(percentiles)$gcm) {
          if (rcm %in% names(data.list[[c(gcm, "uas")]])) {
            ## Select only the values that are within the given season
            i <- data.list[[c(gcm, "uas", rcm, "seasonal.indicies", season)]]

            uas <- data.list[[c(gcm, "uas", rcm, "data")]][i]
            vas <- data.list[[c(gcm, "vas", rcm, "data")]][i]
            pcp <- data.list[[c(gcm, "pr", rcm, "data")]][i]
            temp <- data.list[[c(gcm, "tas", rcm, "data")]][i]
            wind.speed <- data.list[[c(gcm, "wind.speed", rcm, "data")]][i]

            ## Calculate the seasonal thresholds
            temp.thresh <- median(temp)
            pcp.thresh <- quantile(pcp, probs=as.numeric(p))
            wind.speed.thresh <- quantile(wind.speed, probs=as.numeric(p))

            results[p, rcm, season, gcm] <- pinapple.express(uas, vas, pcp, temp, wind.speed, pcp.thresh=pcp.thresh, temp.thresh=temp.thresh, wind.speed.thresh=wind.speed.thresh)["num.events"]
          } else {
            results[p, rcm, season, gcm] <- NA
            warning(paste("Couldn't find rcm ", rcm))
          }
        }
      }
    }
  }
  results
}

apply.blizzard <- function(percentiles, data.list) {

  i <- which(names(dimnames(percentiles)) == "var")
  d <- dim(percentiles)[-i]
  results <- rep(NA, prod(d))
  dim(results) <- d
  dimnames(results) <- dimnames(percentiles)[-i]

  for (p in dimnames(percentiles)$percentiles) {
    for (rcm in dimnames(percentiles)$rcm) {
      for (season in dimnames(percentiles)$season) {
        for (gcm in dimnames(percentiles)$gcm) {
          if (rcm %in% names(data.list[[c(gcm, "pr")]])) {
            ## Select only the values that are within the given season
            i <- data.list[[c(gcm, "pr", rcm, "seasonal.indicies", season)]]

            pcp <- data.list[[c(gcm, "pr", rcm, "data")]][i]
            temp <- data.list[[c(gcm, "tas", rcm, "data")]][i]
            wind.speed <- data.list[[c(gcm, "wind.speed", rcm, "data")]][i]

            ## Calculate the seasonal thresholds
            pcp.thresh <- quantile(pcp, probs=as.numeric(p))
            wind.speed.thresh <- quantile(wind.speed, probs=as.numeric(p))

            results[p, rcm, season, gcm] <- blizzard(pcp, temp, wind.speed, pcp.thresh=pcp.thresh, wind.thresh=wind.speed.thresh)["num.events"]
          } else {
            results[p, rcm, season, gcm] <- NA
            warning(paste("Couldn't find rcm ", rcm))
          }
        }
      }
    }
  }
  results
}

library(RUnit)

CSDI <- function(n, v, years) {
  blocks <- select.blocks.gt.length(d=v, n=n)
  tapply(blocks, years, sum)
}

test.CSDI <- function() {
  f <- CSDI
  cases <- list(list(args=list(n=3, v=c(F, T, T, T, T, F), years=as.factor(rep('2008', 6))),
                     expected=array(4, dimnames=list(c('2008')))),
                list(args=list(n=3, v=c(F, T, T, T, T, F), years=as.factor(c(rep('2008', 2), rep('2009', 4)))),
                     expected=array(c(1, 3), dimnames=list(c('2008', '2009')))),
                list(args=list(n=3, v=c(F, rep(T, 4), rep(F, 2), rep(T, 4)), years=as.factor(c(rep('2008', 6), rep('2009', 5)))),
                     expected=array(c(4, 4), dimnames=list(c('2008', '2009')))),
                list(args=list(n=0, v=c(T, F, T, F, T, F), years=as.factor(rep('2008', 6))),
                     expected=array(3, dimnames=list(c('2008'))))
                )
  for (case in cases) {
    checkEquals(do.call(f, case$args), case$expected)
  }
  checkException(f(n=-1, v=rep(T, 10)))
}
