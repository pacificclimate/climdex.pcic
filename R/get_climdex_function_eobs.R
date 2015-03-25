
#' Returns a list of Climdex variables given constraints
#' Added: monthly only for SPI indices
#' Changed annual only is changed 
#' Added huglin index in tavg although it takes two input variables (tmax and tavg)
#' Changed fre.list to accept also only months and added helper_fun for same reason
get.climdex.variable.list.eobs <- function(source.data.present, time.resolution=c("all", "annual", "monthly"), climdex.vars.subset=NULL) {
  time.res <- match.arg(time.resolution)
  #annual.only <- c("fdETCCDI", "suETCCDI", "idETCCDI", "trETCCDI", "gslETCCDI", "wsdiETCCDI", "csdiETCCDI", "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI", "cddETCCDI", "cwdETCCDI", "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI", "altcddETCCDI", "altcwdETCCDI", "altcsdiETCCDI", "altwsdiETCCDI")
  annual.only <- c("hiETCCDI","gslETCCDI", "wsdiETCCDI", "csdiETCCDI", "cddETCCDI", "cwdETCCDI", "altcddETCCDI", "altcwdETCCDI", "altcsdiETCCDI", "altwsdiETCCDI")
  monthly.only <- c("spi3ETCCDI", "spi6ETCCDI")
  vars.by.src.data.reqd <- list(tmax=c("csuETCCDI", "suETCCDI", "idETCCDI", "txxETCCDI", "txnETCCDI", "tx10pETCCDI", "tx90pETCCDI", "wsdiETCCDI", "altwsdiETCCDI"),
                                tmin=c("cfdETCCDI", "fdETCCDI", "trETCCDI", "tnxETCCDI", "tnnETCCDI", "tn10pETCCDI", "tn90pETCCDI", "csdiETCCDI", "altcsdiETCCDI"),
                                prec=c("spi3ETCCDI", "spi6ETCCDI", "rx1dayETCCDI", "rx5dayETCCDI", "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI",
                                       "cddETCCDI", "cwdETCCDI", "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI", "altcddETCCDI", "altcwdETCCDI"),
                                tavg=c("hd17ETCCDI", "hiETCCDI", "gslETCCDI", "dtrETCCDI") )
  
  if(any(!(source.data.present %in% c("tmin", "tmax", "tavg", "prec"))))
    stop("Invalid variable listed in source.data.present.")
  
  if(all(c("tmax", "tmin") %in% source.data.present) && !("tavg" %in% source.data.present))
    source.data.present <- c(source.data.present, "tavg")
  
  climdex.vars <- unlist(vars.by.src.data.reqd[source.data.present])
  if(!is.null(climdex.vars.subset))
    climdex.vars <- climdex.vars[climdex.vars %in% paste(climdex.vars.subset, "ETCCDI", sep="")]
  
  #freq.lists <- list(c("mon", "yr"), c("yr"))
  freq.lists <- list(c("mon", "yr"), c("yr"),c("mon"))
  
  helper_fun <- function(climdex.vars, annual.only, month.only) {
    if (climdex.vars %in% annual.only) {
      return(paste(climdex.vars, 'yr', sep = '_'))
    } else if (climdex.vars %in% monthly.only) {
      return(paste(climdex.vars, 'mon', sep = '_'))
    } else {
      return(paste(climdex.vars, c('mon', 'yr'), sep = '_'))
    }
  }
  
#   dat <- switch(time.res,
#                 all=unlist(lapply(climdex.vars, function(x) { paste(x, freq.lists[[(x %in% annual.only) + 1]], sep="_") })),
#                 annual=paste(climdex.vars, "yr", sep="_"),
#                 monthly=paste(climdex.vars[!(climdex.vars %in% annual.only)], "mon", sep="_"))
  
  dat <- switch(time.res,
                all=unlist(lapply(climdex.vars, helper_fun, annual.only = annual.only, month.only = month.only)),
                annual=paste(climdex.vars[!(climdex.vars %in% monthly.only)], "yr", sep="_"),
                monthly=paste(climdex.vars[!(climdex.vars %in% annual.only)], "mon", sep="_"))
  
  names(dat) <- NULL
  
  return(dat)
}


#' Returns a list of Climdex functions, with parameters curried in.
#' Create func.names.eca, options.eca and func.eca
#' 
get.climdex.functions.eobs <- function(vars.list, fclimdex.compatible=TRUE) {
  func.names <- c("climdex.fd", "climdex.su", "climdex.id", "climdex.tr", "climdex.gsl",
                  "climdex.fd", "climdex.su", "climdex.id", "climdex.tr",
                  
                  "climdex.txx", "climdex.tnx", "climdex.txn", "climdex.tnn", "climdex.tn10p", "climdex.tx10p", "climdex.tn90p", "climdex.tx90p",
                  "climdex.txx", "climdex.tnx", "climdex.txn", "climdex.tnn", "climdex.tn10p", "climdex.tx10p", "climdex.tn90p", "climdex.tx90p",
                  
                  "climdex.wsdi", "climdex.csdi",
                  
                  "climdex.dtr", "climdex.rx1day", "climdex.rx5day",
                  "climdex.dtr", "climdex.rx1day", "climdex.rx5day",
                  
                  "climdex.sdii", "climdex.r10mm", "climdex.r20mm", "climdex.rnnmm", "climdex.cdd", "climdex.cwd", "climdex.r95ptot", "climdex.r99ptot", "climdex.prcptot",
                  "climdex.sdii", "climdex.r10mm", "climdex.r20mm", "climdex.rnnmm", "climdex.r95ptot", "climdex.r99ptot", "climdex.prcptot",
                  
                  "climdex.cdd", "climdex.cwd", "climdex.csdi", "climdex.wsdi")
  
  func.names.eca <- c("eca.csu", "eca.csu", "eca.cfd", "eca.cfd", "eca.hd17", "eca.hd17", "eca.HI", "eca.SPI3", "eca.SPI6")
  
  el <- list()
  af <- list(freq="annual")
  mf <- list(freq="monthly")
  cwdd.opts <- list(spells.can.span.years=TRUE)
  altcwdd.opts <- list(spells.can.span.years=FALSE)
  wcsdi.opts <- list(spells.can.span.years=FALSE)
  altwcsdi.opts <- list(spells.can.span.years=TRUE)
  rx5day.opts <- list(center.mean.on.last.day=fclimdex.compatible)
  r1mm.opts <- list(threshold=1)
  
  options.climdex <- list( af, af, af, af, el,
                           mf, mf, mf, mf,
                           
                           mf, mf, mf, mf, mf, mf, mf, mf,
                           af, af, af, af, af, af, af, af,
                           
                           wcsdi.opts, wcsdi.opts,
                           
                           mf, mf, c(mf, rx5day.opts),
                           af, af, c(af, rx5day.opts),
                           
                           af, af, af, c(af, r1mm.opts), cwdd.opts, cwdd.opts, af, af, af,
                           mf, mf, mf, c(mf, r1mm.opts), mf, mf, mf,
                           altcwdd.opts, altcwdd.opts, altwcsdi.opts, altwcsdi.opts)
  
  
  options.eca <- list(af,mf, af, mf, af, mf, af, mf, mf)
  
  func.eca <- lapply(1:(length(func.names.eca)), function(n) do.call(functional::Curry, c(list(getFromNamespace(func.names.eca[n], 'ecaclimdex')), options.eca[[n]])))
  
  func.climdex <- lapply(1:(length(func.names)), function(n) do.call(functional::Curry, c(list(getFromNamespace(func.names[n], 'climdex.pcic')), options.climdex[[n]])))
  
  func <- c(func.eca,func.climdex)
  
  names(func) <- c("csuETCCDI_yr", "csuETCCDI_mon", "cfdETCCDI_yr", "cfdETCCDI_mon", "hd17ETCCDI_yr", "hd17ETCCDI_mon", "hiETCCDI_yr", "spi3ETCCDI_mon", "spi6ETCCDI_mon", "fdETCCDI_yr", "suETCCDI_yr", "idETCCDI_yr","trETCCDI_yr", "gslETCCDI_yr",
                   "fdETCCDI_mon", "suETCCDI_mon", "idETCCDI_mon","trETCCDI_mon",
                   
                   "txxETCCDI_mon", "tnxETCCDI_mon", "txnETCCDI_mon", "tnnETCCDI_mon", "tn10pETCCDI_mon", "tx10pETCCDI_mon", "tn90pETCCDI_mon", "tx90pETCCDI_mon",
                   "txxETCCDI_yr",  "tnxETCCDI_yr",  "txnETCCDI_yr",  "tnnETCCDI_yr",  "tn10pETCCDI_yr",  "tx10pETCCDI_yr",  "tn90pETCCDI_yr",  "tx90pETCCDI_yr",
                   
                   "wsdiETCCDI_yr", "csdiETCCDI_yr",
                   
                   "dtrETCCDI_mon", "rx1dayETCCDI_mon", "rx5dayETCCDI_mon",
                   "dtrETCCDI_yr", "rx1dayETCCDI_yr", "rx5dayETCCDI_yr",
                   
                   "sdiiETCCDI_yr", "r10mmETCCDI_yr", "r20mmETCCDI_yr", "r1mmETCCDI_yr", "cddETCCDI_yr", "cwdETCCDI_yr", "r95pETCCDI_yr", "r99pETCCDI_yr", "prcptotETCCDI_yr",
                   "sdiiETCCDI_mon", "r10mmETCCDI_mon", "r20mmETCCDI_mon", "r1mmETCCDI_mon", "r95pETCCDI_mon", "r99pETCCDI_mon", "prcptotETCCDI_mon",
                   
                   "altcddETCCDI_yr", "altcwdETCCDI_yr", "altcsdiETCCDI_yr", "altwsdiETCCDI_yr")
  
  return(func[vars.list])
}


#' Returns metadata for specified Climdex variables
#' Add details for additional indices
get.climdex.variable.metadata.eobs <- function(vars.list, template.filename) {
  all.data <- data.frame(long.name=c("Annual Consecutive Summer Days", "Monthly Consecutive Summer Days",
                                     "Annual Consecutive Frost Days", "Monthly Consecutive Frost Days", 
                                     "Annual Heating Degree days", "Monthly Heating Degree days", "Huglin Index", 
                                     "Standardized Precipitation Index 3mon", 
                                     "Standardized Precipitation Index 6mon",
                                     "Annual Number of Frost Days", "Annual Number of Summer Days", "Annual Number of Icing Days", "Annual Number of Tropical Nights", "Growing Season Length",
                                     "Monthly Number of Frost Days", "Monthly Number of Summer Days", "Monthly Number of Icing Days", "Monthly Number of Tropical Nights",
                                     
                                     "Monthly Maximum of Daily Maximum Temperature", "Monthly Maximum of Daily Minimum Temperature",
                                     "Monthly Minimum of Daily Maximum Temperature", "Monthly Minimum of Daily Minimum Temperature",
                                     
                                     "Percentage of Days when Daily Minimum Temperature is Below the 10th Percentile", "Percentage of Days when Daily Maximum Temperature is Below the 10th Percentile",
                                     "Percentage of Days when Daily Minimum Temperature is Above the 90th Percentile", "Percentage of Days when Daily Maximum Temperature is Above the 90th Percentile",
                                     "Annual Maximum of Daily Maximum Temperature", "Annual Maximum of Daily Minimum Temperature",
                                     "Annual Minimum of Daily Maximum Temperature", "Annual Minimum of Daily Minimum Temperature",
                                     "Percentage of Days when Daily Minimum Temperature is Below the 10th Percentile", "Percentage of Days when Daily Maximum Temperature is Below the 10th Percentile",
                                     "Percentage of Days when Daily Minimum Temperature is Above the 90th Percentile", "Percentage of Days when Daily Maximum Temperature is Above the 90th Percentile",
                                     
                                     "Warm Spell Duration Index", "Cold Spell Duration Index",
                                     
                                      "Mean Diurnal Temperature Range", "Monthly Maximum 1-day Precipitation", "Monthly Maximum Consecutive 5-day Precipitation",
                                     "Mean Diurnal Temperature Range", "Annual Maximum 1-day Precipitation", "Annual Maximum Consecutive 5-day Precipitation",
                                     
                                     "Annual Simple Precipitation Intensity Index", "Annual Count of Days with At Least 10mm of Precipitation",
                                     "Annual Count of Days with At Least 20mm of Precipitation", "Annual Count of Days with At Least 1mm of Precipitation",
                                     "Maximum Number of Consecutive Days with Less Than 1mm of Precipitation", "Maximum Number of Consecutive Days with At Least 1mm of Precipitation",
                                     "Annual Total Precipitation when Daily Precipitation Exceeds the 95th Percentile of Wet Day Precipitation",
                                     "Annual Total Precipitation when Daily Precipitation Exceeds the 99th Percentile of Wet Day Precipitation", "Annual Total Precipitation in Wet Days",
                                     
                                     "Monthly Simple Precipitation Intensity Index", "Monthly Count of Days with At Least 10mm of Precipitation",
                                     "Monthly Count of Days with At Least 20mm of Precipitation", "Monthly Count of Days with At Least 1mm of Precipitation",
                                     "Monthly Total Precipitation when Daily Precipitation Exceeds the 95th Percentile of Wet Day Precipitation",
                                     "Monthly Total Precipitation when Daily Precipitation Exceeds the 99th Percentile of Wet Day Precipitation", "Monthly Total Precipitation in Wet Days",
                                     
                                     "Maximum Number of Consecutive Days Per Year with Less Than 1mm of Precipitation", "Maximum Number of Consecutive Days Per Year with At Least 1mm of Precipitation",
                                     "Cold Spell Duration Index Spanning Years", "Warm Spell Duration Index Spanning Years"),
                       
                         var.name=c("csuETCCDI", "csuETCCDI", "cfdETCCDI", "cfdETCCDI", "hd17ETCCDI", "hd17ETCCDI", "hiETCCDI", "spi3ETCCDI", "spi6ETCCDI",
                                    
                                    "fdETCCDI", "suETCCDI","idETCCDI", "trETCCDI", "gslETCCDI",
                                    "fdETCCDI", "suETCCDI","idETCCDI", "trETCCDI",
                                    
                                    "txxETCCDI", "tnxETCCDI", "txnETCCDI", "tnnETCCDI", "tn10pETCCDI", "tx10pETCCDI", "tn90pETCCDI", "tx90pETCCDI",
                                    "txxETCCDI", "tnxETCCDI", "txnETCCDI", "tnnETCCDI", "tn10pETCCDI", "tx10pETCCDI", "tn90pETCCDI", "tx90pETCCDI",
                                    
                                    "wsdiETCCDI", "csdiETCCDI",
                                    
                                    "dtrETCCDI", "rx1dayETCCDI", "rx5dayETCCDI",
                                    "dtrETCCDI", "rx1dayETCCDI", "rx5dayETCCDI",
                                    
                                    "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI", "cddETCCDI", "cwdETCCDI",  "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI",
                                    "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI", "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI",
                                    
                                    "altcddETCCDI", "altcwdETCCDI", "altcsdiETCCDI", "altwsdiETCCDI"),
                         
                         units=c("days", "days","days","days", "degrees_C", "degrees_C", "degrees_C", "", "", "days", "days", "days","days", "days",
                                 "days", "days", "days","days",
                                 
                                 "degrees_C", "degrees_C", "degrees_C", "degrees_C", "%", "%", "%", "%",
                                 "degrees_C", "degrees_C", "degrees_C", "degrees_C", "%", "%", "%", "%",
                                 
                                 "days", "days",
                                 
                                 "degrees_C", "mm", "mm",
                                 "degrees_C", "mm", "mm",
                                 
                                 "mm d-1", "days", "days", "days", "days", "days", "mm", "mm", "mm",
                                 "mm d-1", "days", "days", "days", "mm", "mm", "mm",
                                 
                                 "days", "days", "days", "days"),
                        
                         annual=c(T, F, T, F, T, F, T, F, F, T, T, T, T, T,
                                  F, F, F, F,
                                  
                                  F, F, F, F, F, F, F, F,
                                  T, T, T, T, T, T, T, T,
                                  
                                  T, T,
                                  
                                  F, F, F,
                                  T, T, T,
                                  
                                  T, T, T, T, T, T, T, T, T,
                                  F, F, F, F, F, F, F,
                                  
                                  T, T, T, T),
                       
                      base.period.attr=c(F, F, F, F, F, F, F, F, F, F, F, F, F, F,
                                            F, F, F, F,
                                            
                                            F, F, F, F, T, T, T, T,
                                            F, F, F, F, T, T, T, T,
                                            
                                            T, T,
                                            
                                            F, F, F,
                                            F, F, F,
                                            
                                            F, F, F, F, F, F, T, T, F,
                                            F, F, F, F, T, T, F,
                                            
                                            F, F, T, T),
                         row.names=c("csuETCCDI_yr", "csuETCCDI_mon", "cfdETCCDI_yr", "cfdETCCDI_mon", "hd17ETCCDI_yr", "hd17ETCCDI_mon", "hiETCCDI_yr", 
                                     "spi3ETCCDI_mon", "spi6ETCCDI_mon",
                                     "fdETCCDI_yr", "suETCCDI_yr", "idETCCDI_yr", "trETCCDI_yr", "gslETCCDI_yr",
                                     "fdETCCDI_mon", "suETCCDI_mon", "idETCCDI_mon", "trETCCDI_mon",
                                     
                                     "txxETCCDI_mon", "tnxETCCDI_mon", "txnETCCDI_mon", "tnnETCCDI_mon", "tn10pETCCDI_mon", "tx10pETCCDI_mon", "tn90pETCCDI_mon", "tx90pETCCDI_mon",
                                     "txxETCCDI_yr", "tnxETCCDI_yr", "txnETCCDI_yr", "tnnETCCDI_yr", "tn10pETCCDI_yr", "tx10pETCCDI_yr", "tn90pETCCDI_yr", "tx90pETCCDI_yr",
                                     
                                     "wsdiETCCDI_yr", "csdiETCCDI_yr",
                                     
                                     "dtrETCCDI_mon", "rx1dayETCCDI_mon", "rx5dayETCCDI_mon",
                                     "dtrETCCDI_yr", "rx1dayETCCDI_yr", "rx5dayETCCDI_yr",
                                     
                                     "sdiiETCCDI_yr", "r10mmETCCDI_yr", "r20mmETCCDI_yr", "r1mmETCCDI_yr", "cddETCCDI_yr", "cwdETCCDI_yr", "r95pETCCDI_yr", "r99pETCCDI_yr", "prcptotETCCDI_yr",
                                     "sdiiETCCDI_mon", "r10mmETCCDI_mon", "r20mmETCCDI_mon", "r1mmETCCDI_mon", "r95pETCCDI_mon", "r99pETCCDI_mon", "prcptotETCCDI_mon",
                                     
                                     "altcddETCCDI_yr", "altcwdETCCDI_yr", "altcsdiETCCDI_yr", "altwsdiETCCDI_yr"),
                         
                         stringsAsFactors=FALSE)
  
  standard.name.lookup <- c(csuETCCDI="consecutive summer days", cfdETCCDI="consecutive frost days", hd17ETCCDI="heating degree days", hiETCCDI="huglin index", 
                            spi3ETCCDI="standardized precipitation index",spi6ETCCDI="standardized precipitation index", 
                            fdETCCDI="number_frost_days", suETCCDI="number_summer_days", idETCCDI="number_icing_days", trETCCDI="number_tropical_nights", gslETCCDI="growing_season_length",
                            txxETCCDI="maximum_daily_maximum_temperature", tnxETCCDI="maximum_daily_minimum_temperature", txnETCCDI="minimum_daily_maximum_temperature", 
                            tnnETCCDI="minimum_daily_minimum_temperature",
                            tn10pETCCDI="percent_days_when_daily_minimum_temperature_below_10p", tx10pETCCDI="percent_days_when_daily_maximum_temperature_below_10p",
                            tn90pETCCDI="percent_days_when_daily_minimum_temperature_above_90p", tx90pETCCDI="percent_days_when_daily_maximum_temperature_above_90p",
                            wsdiETCCDI="warm_spell_duration_index", csdiETCCDI="cold_spell_duration_index", dtrETCCDI="diurnal_temperature_range",
                            altwsdiETCCDI="warm_spell_duration_index", altcsdiETCCDI="cold_spell_duration_index",
                            rx1dayETCCDI="maximum_1day_precipitation", rx5dayETCCDI="maximum_5day_precipitation", sdiiETCCDI="simple_precipitation_intensity_index",
                            r10mmETCCDI="count_days_more_than_10mm_precipitation", r20mmETCCDI="count_days_more_than_20mm_precipitation", r1mmETCCDI="count_days_more_than_1mm_precipitation",
                            cddETCCDI="maximum_number_consecutive_dry_days", cwdETCCDI="maximum_number_consecutive_wet_days",
                            altcddETCCDI="maximum_number_consecutive_dry_days", altcwdETCCDI="maximum_number_consecutive_wet_days",
                            r95pETCCDI="total_precipitation_exceeding_95th_percentile", r99pETCCDI="total_precipitation_exceeding_99th_percentile", prcptotETCCDI="total_wet_day_precipitation")
  
  all.data$standard.name <- standard.name.lookup[all.data$var.name]
  
  all.data$filename <- create.climdex.eobs.filenames(get.split.filename.eobs(template.filename), rownames(all.data))
  return(all.data[vars.list,])
}
