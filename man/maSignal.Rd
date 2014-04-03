\name{maSignal}
\alias{maSignal}
\title{Generate MA plots with or without annotation using LODR estimates}
\usage{
maSignal(exDat, alphaPoint = 0.8, r_mAdjust = T, replicate = T)
}
\arguments{
  \item{exDat}{list, contains input data and stores
  analysis results}

  \item{alphaPoint}{numeric value, for alpha (transparency)
  for plotted points, range is 0 - 1}

  \item{r_mAdjust}{default is TRUE, if FALSE then the r_m
  estimate will not used to offset dashed lines for
  empirical ratios on figure}

  \item{replicate}{default is TRUE, if FALSE then error
  bars will not be produced}
}
\description{
Generate MA plots with or without annotation using LODR
estimates
}

