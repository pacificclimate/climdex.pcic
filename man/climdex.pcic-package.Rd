\name{climdex.pcic}
\alias{climdex.pcic-package}
\docType{package}
\title{climdex.pcic, an implementation of the ETCCDI climate change indices.}
\description{
  This package implements the ETCCDI's 27 core climate change indices
  efficiently in R.
}
\details{
  The calculation of climate extremes are important in many
  disciplines. Annual maximum daily precipitation, annual maximum wind
  speed, and other such extremes are used in many engineering
  applications. However, they are not as useful when speaking about
  climate change.

  The Expert Team on Climate Change Detection and Indicies (ETCCDI)
  has created a set of 27 core indices with the intent of
  capturing the change in the extremes of climate and in selected
  parameters deemed relevant to other disciplines. These model the
  following types of parameters:
  \itemize{
    \item Shifts in the number of days where comparatively extreme
    conditions are observed.
    \item Growing season length.
    \item 10th and 90th percentiles of temperature versus baseline
    conditions.
    \item Lengths of warm, cold, wet, and dry spells.
    \item Counts of days where precipitation exceeds a threshold.
    \item Total precipitation where precipitation exceeds the 95th or
    99th percentile of the baseline.
  }
  
  The \code{climdex.pcic} package provides an implementation of the
  ETCCDI's 27 core climate change indices. It aims to be reasonably high
  performance, to handle non-Gregorian calendar types, to be as correct
  as possible given the definitions of the indices, and to have
  sufficiently readable and concise code as to facilitate easy
  verification by inspection.
}
\references{
  \url{http://cccma.seos.uvic.ca/ETCCDMI/list_27_indices.shtml}
  
  Karl, T.R., N. Nicholls, and A. Ghazi, 1999: CLIVAR/GCOS/WMO
  workshop on indices and indicators for climate extremes: Workshop
  summary. Climatic Change, 42, 3-7.

  Peterson, T.C., and Coauthors: Report on the Activities of the
  Working Group on Climate Change Detection and Related
  Rapporteurs 1998-2001. WMO, Rep. WCDMP-47, WMO-TD 1071, Geneve,
  Switzerland, 143pp.

  Zhang, X., 2005: Avoiding inhomogeneity in percentile-based indices of
  temperature extremes. Journal of Climate 18.11 (2005):1641-.
}
\seealso{
\code{\link{climdexInput.raw}}, \code{\link{climdexInput.csv}}, \code{\link{climdexInput-class}}.
}
\keyword{climate}
\keyword{ts}
