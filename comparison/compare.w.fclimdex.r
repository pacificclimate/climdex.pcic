source('../climdex.r')

station <- '1098D90'

data.dir <- '/home/data/projects/data_cleanup/CDCD_2007/new_data/'
vars <- c(tmax='MAX_TEMP', tmin='MIN_TEMP', prec='ONE_DAY_PRECIPITATION')
fc.data.dir <- '/home/gbuerger/eds/data/PrGeorge/pdd/cdx/'

data.files <- file.path(data.dir, vars, paste(station, '_', vars, '.csv', sep=''))

args <- append(data.files, list(data.column=as.list(vars), base.range=c(1981, 2000)))
clim.in <- do.call(climdexInput, args)

## returns a timeseries object
read.fclimdex <- function(file) {
  data <- read.csv(file, na.strings='-99.9', strip.white=T, sep='', row.names='year')
  ## It's a monthly index
  if ('jan' %in% colnames(data)) {
    year <- as.numeric(rownames(data)[1])
    mon <- 1; stopifnot(colnames(data)[1] == 'jan')
    ts(c(unlist(t(data[,1:12]))), freq=12, start=c(year, mon))
  } else { ## It's an annual index
    year <- as.numeric(rownames(data)[1])
    ts(data, freq=1, start=year)
  }
}

pcic.climdex.to.ts <- function(pc) {
  ## monthly 
  if (all(grepl('[0-9]{4}-[0-9]{2}', names(pc)))) {
    ts(pc, frequency=12, start=as.numeric(strsplit(names(pc)[1], '-', fixed=T)[[1]]))
  } else { # annual
    ts(pc, frequency=1, start=as.numeric(c(names(pc)[1], 1)))
  }
}

compare.pcic.vs.fortran <- function(index) {
  fc.data.file <- file.path(fc.data.dir, toupper(index), paste(station, '.dat', sep=''))

  ## Try to find it if the case is mismatched
  if (! file.exists(fc.data.file)) {
    candidate.files <- list.files(fc.data.dir, paste(station, '.dat', sep=''), recursive=T, full.names=T)
    fc.data.file <- candidate.files[grep(index, candidate.files, ignore.case=T)]
  }

  fc.clim <- read.fclimdex(fc.data.file)
  args <- list(clim.in)

  ## rnnmm is the only index which takes two args
  if (index == 'rnnmm')
    args <- append(args, list(threshold=25))

  pc.clim <- do.call(paste('climdex', index, sep='.'), args)
  pc.clim <- pcic.climdex.to.ts(pc.clim)
  stopifnot(class(pc.clim) == 'ts')

  agree <- fc.clim == pc.clim
  agree.multiplier <- c(NA, 1)[1+ as.numeric(c(!is.na(pc.clim)) & c(agree))]
  agree <- ts(pc.clim * agree.multiplier, start=start(pc.clim), frequency=frequency(pc.clim))
  
  title <- 'Comparison of FClimdex vs. pcic_climdex (R)'
  sub <- paste(toupper(index), station, sep='/')

  pdf(paste(index, '.pdf', sep=''))
  plot(pc.clim, type='b', col='red', main=title, sub=sub, xlab='time', ylab=index)
  lines(fc.clim, type='b')
  points(agree, col='blue')
  legend("top", NULL, c('pcic_climdex', 'FClimdex'), lty=c(1, 1, NA), pch=c(NA, NA, 1), col=c('red', 'black', 'blue'))
  dev.off()
}

for (index in all.indicies) {
  rv <- try(compare.pcic.vs.fortran(index))
  ## Cleanup
  if (inherits(rv, 'try-error')) {
    dev.off()
    print(paste("Failed to compare index", index))
  }
}
