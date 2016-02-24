fname = '~/Downloads/TX_STAID000162.txt'

find_start_of_data_index = function(fname) {
  matched_indices = which(grepl('^SOUID', gsub(" ", "", readLines(fname, n = 50))))
  if (length(matched_indices) > 1) stop('ECA fileformat error: cannot determine start of data, multiple header lines')
  if (length(matched_indices) == 0) stop('ECA fileformat error: cannot find start of data: cannot find header line')
  return(matched_indices - 1)
}


tg <- eca.input('~/Documents/deBilt_tavg/TG_STAID000162.txt', 'TG', 'DATE')
tx <- eca.input('~/Documents/deBilt_tmax/TX_STAID000162.txt', 'TX', 'DATE')
tn <- eca.input('~/Documents/deBilt_tmin/TN_STAID000162.txt', 'TN', 'DATE')
prec <- eca.input('~/Documents/deBilt_prec/RR_STAID000162.txt', 'RR', 'DATE')
fx <- eca.input('~/Documents/deBilt_windgust/FX_STAID000162.txt', 'FX', 'DATE')

ci <- climdexInput.raw(tmax= tx$TX, tmin=tn$TN, tavg=tg$TG, prec=prec$RR, wind_gust = fx$FX, 
                       tmax.dates = tx$DATE, tmin.dates = tn$DATE, tavg.dates = tg$DATE, prec.dates = prec$DATE, wind_gust.dates = fx$DATE,
                       base.range=c(1961, 1991))


ci <- climdexInput.raw(tmax=tx$TX, tavg=tg$TG, prec=prec$RR, wind_gust = fx$FX, 
                       tmax.dates = tx$DATE, tavg.dates = tg$DATE, prec.dates = prec$DATE, wind_gust.dates = fx$DATE)



ci <- climdexInput.raw(tavg=tg$TG, prec=prec$RR, wind_gust = fx$FX, 
                        tavg.dates = tg$DATE, prec.dates = prec$DATE, wind_gust.dates = fx$DATE)


test_max_fx <- climdex.wind_maxgust(ci, "halfyear")

test_rx1day <- climdex.rx1day(ci, "halfyear")

test_cd <- climdex.cd(ci, 'annual')


rm(list = ls())
ci <- climdexInput.raw(prec=prec$RR, prec.dates = prec$DATE, base.range=c(1961,1990))


