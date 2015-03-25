#' STANDARDIZE PRECIPITATION INDEX
#' based on the formula in the spi package 
#' SPI is only has monthly outputs (see changes in get.climdex.variable.list.eobs)
#' make monthly sums
monthly_sums_spi <- function(temp,date.factor) {
  stopifnot(is.numeric(temp) && is.factor(date.factor))
  return(tapply(temp,date.factor,sum,na.rm=TRUE))
}
#' calcualte eca monthly sums
eca_sums_mon <- function(ci,freq=c("monthly")){
  stopifnot(!is.null(ci@data$prec))
  return(monthly_sums_spi(ci@data$prec,ci@date.factors$monthly) * ci@namasks$monthly$prec)
}
## Makes multi-month averages depending on k (here k=3,6). For spi3 k=3
# precipitation is a vector of monthly precipitation sums
# returns monthly precipation averaged over current month and prior k-1 months
getPrecOnTimescale <- function(precipitation,k){
  Nt <- length(precipitation)
  prec.k <- as.vector(sapply(seq(from=1, to=Nt),function(t) {tm <- max(t-k+1,1); sum(as.vector(precipitation[tm:t]))}))
  return(prec.k)
}

#'In ECA&D if 85% of data are present then we consider them as valid. 
#'This was not used in SPI or elsewhere (yet)
# dataLeadsToNAindex <- function(dat) {
#   percentage_should_be_present = 0.85
#   number_of_nan = sum(is.na(dat))
#   total_length = length(dat)
#   return((number_of_nan / total_length) < percentage_should_be_present)
# }

#' empirical Standard Precipitation Index 3 (k in getPrecOnTimescale is 3)
eca.SPI3 <- function(ci,freq=c("monthly")){
  dat_mon <- eca_sums_mon(ci, freq=c("monthly"))
  dat <- getPrecOnTimescale(dat_mon,3)
  #if(dataLeadsToNAindex(dat)) return(NA)
  if(all(is.na(dat))) return(NA)
  
  fit.cdf <- ecdf(dat)
  cdfs <- fit.cdf(dat)
  spi.t <- qnorm(cdfs)
  spi.sym <- spi.t
  # drop Inf
  spi.sym[which(spi.t == Inf)] <- NA
  spi.sym[which(spi.t == -Inf)] <-NA
  return(spi.sym)
}

#' empirical Standard Precipitation Index 6 (k in getPrecOnTimescale is 6)
eca.SPI6 <- function(ci,freq=c("monthly")){
  dat_mon <- eca_sums_mon(ci, freq=c("monthly"))
  dat <- getPrecOnTimescale(dat_mon,6)
  #if(dataLeadsToNAindex(dat)) return(NA)
  if(all(is.na(dat))) return(NA)
  
  fit.cdf <- ecdf(dat)
  cdfs <- fit.cdf(dat)
  spi.t <- qnorm(cdfs)
  spi.sym <- spi.t
  # drop Inf (in the spi package Inf is set to 2.2 and -Inf to -2.2)
  spi.sym[which(spi.t == Inf)] <- NA
  spi.sym[which(spi.t == -Inf)] <-NA
  return(spi.sym)
}


#' CONSECUTIVE SUMMER DAYS
#' maximum number of concsecutive summer days (tx > 25C)
#' spell.can.span.years for us is FALSE
eca.csu <- function(ci,freq=c("monthly","annual"),spells.can.span.years=FALSE) {
  stopifnot(!is.null(ci@data$tmax));
  return(climdex.pcic:::spell.length.max(ci@data$tmax, ci@date.factors[[match.arg(freq)]], 25, ">", spells.can.span.years) * ci@namasks[[match.arg(freq)]]$tmax)
}

#' CONSECUTIVE FROST DAYS
#' maximum number of concsecutive frost days (tn < 0C)
#' spell.can.span.years for us is FALSE
eca.cfd <- function(ci,freq=c("monthly","annual"),spells.can.span.years=FALSE) { 
  stopifnot(!is.null(ci@data$tmin));
  return(climdex.pcic:::spell.length.max(ci@data$tmin, ci@date.factors[[match.arg(freq)]], 0, "<", spells.can.span.years) * ci@namasks[[match.arg(freq)]]$tmin) }

#' CONSECUTIVE WET DAYS
#' We would This is a trial to have cwd for months also (not sure about this)
#'  spells.can.span.years=FALSE
climdex.cwd.eobs <- function(ci, freq=c("monthly","annual"),spells.can.span.years=FALSE) {
  stopifnot(!is.null(ci@data$prec)); 
  return(climdex.pcic:::spell.length.max(ci@data$prec, ci@date.factors[[match.arg(freq)]], 1, ">=", spells.can.span.years) * ci@namasks[[match.arg(freq)]]$prec)
}

#' HEATING DEGREE DAYS 
#' sum of (17C-TG) (C)
#' 
eca.hd17 <- function(ci, freq=c("monthly","annual")) { 
  stopifnot(!is.null(ci@data$tavg)); 
  return(tapply((17 -  ci@data$tavg), ci@date.factors[[match.arg(freq)]], sum)* ci@namasks[[match.arg(freq)]]$tavg)
}

######################################################################################################################
#' HUGLIN INDEX
#' This function is not ready yet. It is uses coeeficient based on the latitude
#' For this index I had to curry the cdx.funcs to be able to include the subset. Later I realised 
#' I need also to include the latitude. This would be used in compute.indices.for.stripe together with get.lat 
#' to retrieve subset & latitude 
#' I didn't proceed with finisheing the eca.HI function. I thought to ask you first if its possible 
#' to adapt compute.indices.for.stripe so it can include the currying and the latitude. Or if you had a better idea on this 
#' please let me know. 

#' Function to Curry a cxd.funcs for subset
#' used only for Huglin Index
curry_in_subset_for_huglin <- function(cdx.funcs, subset){
  cdx.names = names(cdx.funcs)
  cdx.funcs <- lapply(cdx.names, function(function_name) {
    if(grepl('^hi', function_name)) {
      f = cdx.funcs[[function_name]]
      return(functional::Curry(f, subset = subset[['Y']]))
    } else {
      return(f)
    }
  })
  names(cdx.funcs) = cdx.names
  return(cdx.funcs)
}

#' Function to retrieve the latitude 
#'  we have different variable names for temp & precip; that's why I've changed the variable.name.map
get.lat <- function(open_file_list, v.f.idx, variable.name.map=c(tmax="tx", tmin="tn", prec="rr", tavg="tg")) {
  var.name <- variable.name.map[[names(v.f.idx)[1]]]
  y.dim <- ncdf4.helpers::nc.get.dim.for.axis(open_file_list[[1]], var.name, "Y")
  return(y.dim$vals)
}

#' Huglin index 
#' It is annual only and valid for months (April-Sep)
eca.HI <- function(ci,freq=c("annual"),subset=-1){
  
  tempavg <- ci@data$tavg
  tempmax <- ci@data$tmax  
  
  #browser()
  month.series <- climdex.pcic:::get.months(ci@dates)
  year.series <- climdex.pcic:::get.years(ci@dates)
  valid.months <- month.series >=4 & month.series <=9
  #browser()
  hi_coef <-  if (subset <=40) {hi_coeff <- 1 
  }else if(subset >40 & subset <42) {hi_coef <- 1.02
  }else if(subset >42 & subset <44) {hi_coef <- 1.03
  }else if(subset >44 & subset <48) {hi_coef <- 1.04
  }else if(subset >46 & subset <48) {hi_coef <- 1.05
  }else if(subset >48 & subset <50) {hi_coef <- 1.06 
  }else if(subset >=50){hi_coef <- 1}
  valid.sel<- year.series[valid.months]
  tempdata <- ((((tempavg -10) + (tempmax -10)) /2) * hi_coef)
  dat_final <- tempdata[valid.months]
  return(tapply(dat_final,valid.sel,sum,na.rm=TRUE))
  
}
