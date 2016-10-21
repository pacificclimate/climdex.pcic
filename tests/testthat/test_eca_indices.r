#' Description
#' This testthat tests the new indices added in climdex.pcic 
#' 1. Test the eca.input function which can be used for station files from ECA&D
#' 2. Test the climdexInput.raw for the added parameters, added frequencies (halfyear, seasons), and new quantiles
#' 3. Test new indices for precipitation and temperature
#' 4. test new indices for new parameters: cloud, sun, snow, wind.
#' 5. We don't test yet the sun_rel and the snow_new (related indices)

rm(list = ls())
library(climdex.pcic)
library(testthat)

####################################
### Precipitation indices
####################################
context("Precipitation indices")

test_that("Precipitation indices annual", {
  prec <- eca.input('~/Documents/deBilt_prec/RR_STAID000162.txt', 'RR', 'DATE')
  
  expect_equal_to_reference(prec, "./tests/testthat/outputTests/rr_deBilt.rds")
  
  ## Here we test the additional frequencies (halfyear & seasons)
  ci <- climdexInput.raw(prec=prec$RR, prec.dates = prec$DATE, base.range=c(1961,1990))
  expect_equal_to_reference(ci, "./tests/testthat/outputTests/ci_rr_deBilt.rds")
  
  ## Test new indices for precipitation
  ci_spi3   <- climdex.spi3(ci, freq=c("monthly"), scale=3)
  ci_spi6   <- climdex.spi6(ci, freq=c("monthly"), scale=6)
  ci_r99p <- climdex.r99p(ci, freq=c("monthly"))
  ci_r99ptot <- climdex.r99ptot(ci, freq=c("monthly"))
  
  expect_equal_to_reference(ci_spi3, "./tests/testthat/outputTests/spi3_deBilt.rds")
  expect_equal_to_reference(ci_spi6, "./tests/testthat/outputTests/spi6_deBilt.rds")
  expect_equal_to_reference(ci_r99p, "./tests/testthat/outputTests/r99p_deBilt.rds")
  expect_equal_to_reference(ci_r99ptot, "./tests/testthat/outputTests/r99ptot_deBilt.rds")
  

})


####################################
### Temperature indices
####################################

context("tmax indices")

test_that("tmax indices annual", {
  tg <- eca.input('~/Documents/deBilt_tavg/TG_STAID000162.txt', 'TG', 'DATE')
  tx <- eca.input('~/Documents/deBilt_tmax/TX_STAID000162.txt', 'TX', 'DATE')
  tn <- eca.input('~/Documents/deBilt_tmin/TN_STAID000162.txt', 'TN', 'DATE')
  prec <- eca.input('~/Documents/deBilt_prec/RR_STAID000162.txt', 'RR', 'DATE')
  
  expect_equal_to_reference(tg, "./tests/testthat/outputTests/tg_deBilt.rds")
  expect_equal_to_reference(tx, "./tests/testthat/outputTests/tx_deBilt.rds")
  expect_equal_to_reference(tn, "./tests/testthat/outputTests/tn_deBilt.rds")
  expect_equal_to_reference(prec, "./tests/testthat/outputTests/rr_deBilt.rds")
  
  
  ## Here we test the additional frequencies (halfyear & seasons) and the additional quantiles for temperature (q25, q75 )
  ci_temp <- climdexInput.raw(tmax= tx$TX, tmin=tn$TN, tavg=tg$TG, tmax.dates = tx$DATE, tmin.dates = tn$DATE, 
                         prec=prec$RR, prec.dates = prec$DATE, 
                         tavg.dates = tg$DATE,base.range=c(1961, 1991))
  
  expect_equal_to_reference(ci_temp, "./tests/testthat/outputTests/ci_temp_deBilt.rds")
  
  ## tmax
  ci_csu <- climdex.csu(ci_temp, freq=c("annual"))
  ci_txndaymin <- climdex.txndaymin(ci_temp, freq=c("annual"))
  ci_txndaymax <- climdex.txndaymin(ci_temp, freq=c("annual"))
  
  expect_equal_to_reference(ci_csu, "./tests/testthat/outputTests/csu_deBilt.rds")
  expect_equal_to_reference(ci_txndaymin, "./tests/testthat/outputTests/txndaymin_deBilt.rds")
  expect_equal_to_reference(ci_txndaymax, "./tests/testthat/outputTests/txndaymax_deBilt.rds")
})

context("tmin indices")

test_that("tmin indices annual & monthly", {
  tn <- eca.input('~/Documents/deBilt_tmin/TN_STAID000162.txt', 'TN', 'DATE')

  ci_temp <- climdexInput.raw(tmin=tn$TN, tmin.dates = tn$DATE,base.range=c(1961, 1991))
  
  ci_cfd   <- climdex.cfd(ci_temp, freq=c("monthly"))
  ci_tnndaymin <- climdex.tnndaymin(ci_temp, freq=c("monthly"))
  ci_tnndaymax <- climdex.tnndaymax(ci_temp, freq=c("annual"))
  
  expect_equal_to_reference(ci_cfd, "./tests/testthat/outputTests/cfd_deBilt.rds")
  expect_equal_to_reference(ci_tnndaymin, "./tests/testthat/outputTests/tnndaymin_deBilt.rds")
  expect_equal_to_reference(ci_tnndaymax, "./tests/testthat/outputTests/tnndaymax_deBilt.rds")
})

context("tavg indices")

test_that("tavg indices annual & monthly", {
  tn <- eca.input('~/Documents/deBilt_tmin/TN_STAID000162.txt', 'TN', 'DATE')
  tg <- eca.input('~/Documents/deBilt_tavg/TG_STAID000162.txt', 'TG', 'DATE')
  prec <- eca.input('~/Documents/deBilt_prec/RR_STAID000162.txt', 'RR', 'DATE')
  
  ci_temp <- climdexInput.raw(tmin = tn$TN, tavg = tg$TG, prec = prec$RR, 
                              tmin.dates = tn$DATE, prec.dates = prec$DATE, tavg.dates = tg$DATE,
                              base.range=c(1961, 1991))
  
  ##tavg
  ci_hd17        <- climdex.hd17(ci_temp, freq=c("annual"))
  ci_tmndaymin   <- climdex.tmndaymin(ci_temp, freq=c("monthly"))
  ci_tmndaymax   <- climdex.tmndaymax(ci_temp, freq=c("monthly")) 
  ci_cd        <- climdex.cd(ci_temp, freq=c("annual"))
  ci_cw        <- climdex.cw(ci_temp, freq=c("monthly"))
  ci_wd        <- climdex.wd(ci_temp, freq=c("monthly"))
  ci_ww        <- climdex.ww(ci_temp, freq=c("monthly"))
  
  expect_equal_to_reference(ci_hd17, "./tests/testthat/outputTests/hd17_deBilt.rds")
  expect_equal_to_reference(ci_tmndaymin, "./tests/testthat/outputTests/tmndaymin_deBilt.rds")
  expect_equal_to_reference(ci_tmndaymax, "./tests/testthat/outputTests/tmndaymax_deBilt.rds")
  expect_equal_to_reference(ci_cd, "./tests/testthat/outputTests/cd_deBilt.rds")
  expect_equal_to_reference(ci_cw, "./tests/testthat/outputTests/cw_deBilt.rds")
  expect_equal_to_reference(ci_wd, "./tests/testthat/outputTests/wd_deBilt.rds")
  expect_equal_to_reference(ci_ww, "./tests/testthat/outputTests/ww_deBilt.rds")
  
})

####################################
### Wind indices
####################################

context("Wind speed")

test_that("Wind speed annual & monthly", {
  fg <- eca.input('~/Documents/deBilt_windspeed/FG_STAID000162.txt', 'FG', 'DATE')
  fx <- eca.input('~/Documents/deBilt_windgust/FX_STAID000162.txt', 'FX', 'DATE')
  dd <- eca.input('~/Documents/deBilt_winddir/DD_STAID000162.txt', 'DD', 'DATE')
  
  expect_equal_to_reference(fg, "./tests/testthat/outputTests/fg_deBilt.rds")
  expect_equal_to_reference(fx, "./tests/testthat/outputTests/fx_deBilt.rds")
  expect_equal_to_reference(dd, "./tests/testthat/outputTests/dd_deBilt.rds")

  
  ## Here we test the additional frequencies (halfyear & seasons) and the additional quantiles for temperature (q25, q75 )
  ci_wind <- climdexInput.raw(wind= fg$FG, wind_gust=fx$FX, wind_dir=dd$DD, wind.dates = fg$DATE, wind_gust.dates = fx$DATE, 
                              wind_dir.dates = dd$DATE, base.range=c(1961, 1991))
  
  expect_equal_to_reference(ci_wind, "./tests/testthat/outputTests/ci_wind_deBilt.rds")
  
  ## wind speed
  ci_fg   <- climdex.fg(ci_wind, freq=c("monthly"))
  ci_fgcalm <- climdex.fgcalm(ci_wind, freq=c("monthly"))
  ci_fg6bft <- climdex.fg6bft(ci_wind, freq=c("annual"))
  
  expect_equal_to_reference(ci_fg, "./tests/testthat/outputTests/fg_ci_deBilt.rds")
  expect_equal_to_reference(ci_fgcalm, "./tests/testthat/outputTests/fgcalm_deBilt.rds")
  expect_equal_to_reference(ci_fg6bft, "./tests/testthat/outputTests/fg6bft_deBilt.rds")
  
})
# 
# context("Wind gust")
# 
# test_that("Wind gust monthly", {
#   fx <- eca.input('~/Documents/deBilt_windgust/FX_STAID000162.txt', 'FX', 'DATE')
# 
#   ci_wind <- climdexInput.raw(wind_gust=fx$FX, wind_gust.dates = fx$DATE, base.range=c(1961, 1991))
#   
#   ## wind gust
#   ci_fxstorm   <- climdex.fxstorm(ci_wind, freq=c("monthly"))
#   ci_fxx <- climdex.fxx(ci_wind, freq=c("monthly"))
# 
#   expect_equal_to_reference(ci_fxstorm, "./tests/testthat/outputTests/fxstorm_deBilt.rds")
#   expect_equal_to_reference(ci_fxx, "./tests/testthat/outputTests/fxx_deBilt.rds")
# })
# 
# context("Wind direction")
# 
# test_that("Wind direction annual and monthly", {
#   dd <- eca.input('~/Documents/deBilt_winddir/DD_STAID000162.txt', 'DD', 'DATE')
#   
#   ci_wind <- climdexInput.raw(wind_dir=dd$DD, wind_dir.dates = dd$DATE, base.range=c(1961, 1991))
#   
#   ## wind direction
#   ci_ddnorth   <- climdex.ddnorth(ci_wind, freq=c("monthly"))
#   ci_ddeast <- climdex.ddeast(ci_wind, freq=c("monthly"))
#   ci_ddsouth <- climdex.ddsouth(ci_wind, freq=c("annual"))
#   ci_ddwest <- climdex.ddwest(ci_wind, freq=c("annual"))
#   
#   expect_equal_to_reference(ci_ddnorth, "./tests/testthat/outputTests/ddnorth_deBilt.rds")
#   expect_equal_to_reference(ci_ddeast, "./tests/testthat/outputTests/ddeast_deBilt.rds")
#   expect_equal_to_reference(ci_ddsouth, "./tests/testthat/outputTests/ddsouth_deBilt.rds")
#   expect_equal_to_reference(ci_ddwest, "./tests/testthat/outputTests/ddwest_deBilt.rds")
# 
# })

####################################
### Snow indices
####################################
# 
# context("Snow depth")
# 
# test_that("Snow annual & monthly", {
#   snowD <- eca.input('~/Documents/Lugano_snowdepth/SD_STAID000242.txt', 'SD', 'DATE')
# 
#   expect_equal_to_reference(snowD, "./tests/testthat/outputTests/snowD_deBilt.rds")
#   
#   ci_snow <- climdexInput.raw(snow = snowD$SD, snow.dates = snowD$DATE, base.range=c(1961, 1991))
#   
#   expect_equal_to_reference(ci_snow, "./tests/testthat/outputTests/ci_snow_Lugano.rds")
#   
#   ## snow depth
#   ci_sdd   <- climdex.sdd(ci_snow, freq=c("monthly"))
#   ci_sdx <- climdex.sdx(ci_snow, freq=c("monthly"))
#   ci_sd <- climdex.sd(ci_snow, freq=c("annual"))
#   
#   expect_equal_to_reference(ci_sdd, "./tests/testthat/outputTests/sdd_deBilt.rds")
#   expect_equal_to_reference(ci_sdx, "./tests/testthat/outputTests/sdx_deBilt.rds")
#   expect_equal_to_reference(ci_sd, "./tests/testthat/outputTests/sd_deBilt.rds")
#   
# })

####################################
### Cloud indices
####################################

context("Cloud indices")

test_that("Cloud indices annual & monthly", {
  cloud <- eca.input('~/Documents/deBilt_cloud/CC_STAID000162.txt', 'CC', 'DATE')
  
  expect_equal_to_reference(cloud, "./tests/testthat/outputTests/cloud_deBilt.rds")
  
  ci_cloud <- climdexInput.raw(cloud = cloud$CC, cloud.dates = cloud$DATE, base.range=c(1961, 1991))
  
  expect_equal_to_reference(ci_cloud, "./tests/testthat/outputTests/ci_cloud_deBilt.rds")
  
  ## snow depth
  ci_cc   <- climdex.cc(ci_cloud, freq=c("monthly"))
  ci_cc6 <- climdex.cc6(ci_cloud, freq=c("monthly"))
  ci_cc2 <- climdex.cc2(ci_cloud, freq=c("annual"))
  
  expect_equal_to_reference(ci_cc, "./tests/testthat/outputTests/cc_deBilt.rds")
  expect_equal_to_reference(ci_cc6, "./tests/testthat/outputTests/cc6_deBilt.rds")
  expect_equal_to_reference(ci_cc2, "./tests/testthat/outputTests/cc2_deBilt.rds")
  
})

####################################
### Sun indices
####################################

context("Sun indices")

test_that("Sun index monthly", {
  sunshine <- eca.input('~/Documents/deBilt_sunshine/SS_STAID000162.txt', 'SS', 'DATE')
  
  expect_equal_to_reference(sunshine, "./tests/testthat/outputTests/sunshine_deBilt.rds")
  
  ci_sun <- climdexInput.raw(sun = sunshine$SS, sun.dates = sunshine$DATE, base.range=c(1961, 1991))
  
  expect_equal_to_reference(ci_sun, "./tests/testthat/outputTests/ci_sun_deBilt.rds")
  
  ## snow depth
  ci_ss  <- climdex.ss(ci_sun, freq=c("monthly"))

  expect_equal_to_reference(ci_ss, "./tests/testthat/outputTests/ss_deBilt.rds")

})


