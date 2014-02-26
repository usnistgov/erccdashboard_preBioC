\name{estLODR}
\alias{estLODR}
\title{Estimate Limit of Detection of Ratios (LODR)}
\usage{
estLODR(expDat, kind = "ERCC", prob = 0.9)
}
\arguments{
  \item{expDat}{list, contains input data and stores
  analysis results}

  \item{kind}{"ERCC" or "Sim"}

  \item{prob}{probability, ranging from 0 - 1, default is
  0.9}
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

