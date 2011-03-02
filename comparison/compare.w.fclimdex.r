## compare.w.fclimdex.r
## James Hiebert <hiebert@uvic.ca>
##
## Description: Creates comparison plots of pcic_climdex vs. FClimdex by running pcic_climdex on
## station data and comparing that to an output directory from FClimdex
## Accepted arguments are: station, data.dir, vars, fc.data.dir, and base.range
## See the bottom of the file for default arguments
##
## You may call this script from the shell such as this:
## R CMD BATCH '--args station="1098D90" base.range=c(1971, 2000)' compare.w.f.climdex.r
##
## TOOD: plot.dir argument

source('../climdex.r')

## returns a timeseries object
read.fclimdex <- function(file) {
  data <- read.csv(file, na.strings=c('-99.9', '-99.90'), strip.white=T, sep='', row.names='year')
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

## returns a timeseries object
pcic.climdex.to.ts <- function(pc) {
  ## monthly 
  if (all(grepl('[0-9]{4}-[0-9]{2}', names(pc)))) {
    rv <- ts(matrix(pc), frequency=12, start=as.numeric(strsplit(names(pc)[1], '-', fixed=T)[[1]]))
  } else { # annual
    rv <- ts(matrix(pc), frequency=1, start=as.numeric(c(names(pc)[1], 1)))
  }
}

get.climdex.agreement <- function(pc.clim, fc.clim) {
  agree <- fc.clim == pc.clim
  agree.multiplier <- c(NA, 1)[1+ as.numeric(c(!is.na(pc.clim)) & c(agree))]
  ts(pc.clim * agree.multiplier, start=start(pc.clim), frequency=frequency(pc.clim))
}

## Returns a ts object with T values for when one model is NA but the other is not
get.climdex.na.disagreement <- function(one, other) {
  one@.Data <- is.na(one@.Data)
  other@.Data <- is.na(other@.Data)
  rv <- xor(one, other)
  replace(rv, rv == F, NA)
}

plot.double.climdex <- function(pc.clim, fc.clim, ...) {
  agree <- get.climdex.agreement(pc.clim, fc.clim)

  ylim <- range(c(pc.clim, fc.clim), na.rm=T)
  plot(pc.clim, type='b', col='red', ylim=ylim, ...)
  lines(fc.clim, type='b')
  points(agree, col='blue')
  legend("top", NULL, c('pcic_climdex', 'FClimdex', 'agreement'), lty=c(1, 1, NA), pch=c(1, 1, 1), col=c('red', 'black', 'blue'))
}

plot.diff.climdex <- function(pc.clim, fc.clim, ...) {
  agree <- get.climdex.agreement(pc.clim, fc.clim)
  agree <- replace(agree, !is.na(agree), 0)
  ylim <- range(pc.clim - fc.clim, 0, na.rm=T)
  plot(pc.clim - fc.clim, type='b', col='red', ylim=ylim, ...)

  na.diffs <- get.climdex.na.disagreement(pc.clim, fc.clim)
  points(replace(na.diffs, na.diffs==T, 0), col='red')

  points(agree, col='blue')
  legend("top", NULL, c("pcic_climdex - Fclimdex", "Agreement"), lty=c(1, NA), pch=c(1, 1), col=c('red', 'blue'))
}

compare.pcic.vs.fortran <- function(index, clim.in, station, fc.data.dir) {
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

  title <- 'Comparison of FClimdex vs. pcic_climdex (R)'
  sub <- paste(toupper(index), station, sep='/')

  pdf(paste(index, '.pdf', sep=''))
  plot.double.climdex(pc.clim, fc.clim, main=title, sub=sub, xlab='time', ylab=index)

  title <- 'Difference between pcic_climdex (R) and FClimdex'
  plot.diff.climdex(pc.clim, fc.clim, main=title, sub=sub, xlab='time', ylab=index)
  dev.off()
}

run.it <- function(station, data.dir, vars, fc.data.dir, base.range) {

  data.files <- file.path(data.dir, vars, paste(station, '_', vars, '.csv', sep=''))

  args <- append(data.files, list(data.column=as.list(vars), base.range=base.range))
  clim.in <- do.call(climdexInput, args)

  for (index in all.indicies) {
    rv <- try(compare.pcic.vs.fortran(index, clim.in, station, fc.data.dir))
    ## Cleanup
    if (inherits(rv, 'try-error')) {
      if (dev.cur() != 1) dev.off()
      print(paste("Failed to compare index", index))
    }
  }
}

## Default arguments to the script
default.args <- list(station =     '1098D90',
                     data.dir =    '/home/data/projects/data_cleanup/CDCD_2007/new_data/',
                     vars =        c(tmax='MAX_TEMP', tmin='MIN_TEMP', prec='ONE_DAY_PRECIPITATION'),
                     fc.data.dir = '/home/gbuerger/eds/data/PrGeorge/pdd/cdx/',
                     base.range = c(1981, 2000)
                     )

## Assign given arguments to the current environment
args <- commandArgs(trailingOnly=T)
for (a in args) {
  print(a)
  eval(parse(text=a))
}

## Fill in (from defaults) arguments which were not provided
provided.args <- ls(pattern=paste("(", paste(names(default.args), collapse="|"), ")", sep=''))
arg.list <- lapply(names(default.args),
                   function(n) {if (n %in% provided.args) get(n, parent.frame(2))
                                else assign(n, default.args[[n]], parent.frame(2))}
                   )
names(arg.list) <- names(default.args)

## Invoke the script
run.it(station, data.dir, vars, fc.data.dir, base.range)
