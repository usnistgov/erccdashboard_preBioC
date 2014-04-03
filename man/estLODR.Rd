\name{estLODR}
\alias{estLODR}
\title{Estimate Limit of Detection of Ratios (LODR)}
\usage{
estLODR(exDat, kind = "ERCC", prob = 0.9, yrange = c(1e-15, 10))
}
\arguments{
  \item{exDat}{list, contains input data and stores
  analysis results}

  \item{kind}{"ERCC" or "Sim"}

  \item{prob}{probability, ranging from 0 - 1, default is
  0.9}

  \item{yrange}{y-axis limits, default is c(1e-15, 1e1)}
}
\description{
Estimate Limit of Detection of Ratios (LODR)
}
\details{
This is the function to estimate a limit of detection of
ratios (LODR) for a a chosen probability and threshold
p-value for the fold changes in the ERCC control ratio
mixtures.
}

