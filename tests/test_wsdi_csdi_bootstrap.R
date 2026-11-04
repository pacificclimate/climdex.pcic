library(climdex.pcic)
library(RUnit)

## WSDI and CSDI use the Zhang et al. (2005) in-base bootstrap thresholds for
## years inside the base period, as the TX90p-type indices do.
climdex.pcic.test.wsdi.csdi.bootstrap <- function() {
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, NULL,
                         tmax.dates, tmin.dates, NULL, base.range=c(1981, 1990))

  for(index in c("wsdi", "csdi")) {
    fun <- match.fun(paste("climdex", index, sep="."))
    for(span in c(FALSE, TRUE)) {
      old <- fun(ci, spells.can.span.years=span, bootstrap=FALSE)
      new <- fun(ci, spells.can.span.years=span)
      years <- as.numeric(names(new))
      in.base <- years >= 1981 & years <= 1990

      ## Years outside the base period are unchanged.
      checkEquals(old[!in.base], new[!in.base], paste(index, "changed outside the base period"))
      ## Years inside the base period use the bootstrap thresholds.
      checkTrue(any(abs(old[in.base] - new[in.base]) > 0, na.rm=TRUE), paste(index, "unchanged inside the base period"))
      checkTrue(all(new[in.base] >= 0, na.rm=TRUE))
    }
  }

  ## Check one base year directly against the definition: the mean, over the
  ## (n - 1) bootstrap threshold sets for that year, of the number of days in
  ## spells of at least 6 days.
  q90 <- ci@quantiles$tmax$inbase$q90
  yr <- as.numeric(format(ci@dates, "%Y"))
  idx <- which(yr == 1985)
  expected <- mean(sapply(1:dim(q90)[3], function(k) {
    sum(select.blocks.gt.length(ci@data$tmax[idx] > q90[ci@jdays[idx], 1985 - 1981 + 1, k], 5))
  }))
  checkEqualsNumeric(expected, climdex.wsdi(ci)["1985"])
}
