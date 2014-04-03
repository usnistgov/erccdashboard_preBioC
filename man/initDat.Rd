\name{initDat}
\alias{initDat}
\title{Initialize the exDat list}
\usage{
initDat(datType = NULL, exTable = NULL, repNormFactor = NULL,
  filenameRoot = NULL, sample1Name = NULL, sample2Name = NULL,
  erccmix = "RatioPair", erccdilution = 1, spikeVol = 1,
  totalRNAmass = 1, choseFDR = 0.05, ratioLim = c(-4, 4),
  signalLim = c(-14, 14), userMixFile = NULL)
}
\arguments{
  \item{datType}{type is "count" or "array", unnormalized
  data is expected (normalized data may be accepted in
  future version of the package). Default is "count"
  (integer count data),"array" is unnormalized fluorescent
  intensities from microarray fluorescent intensities (not
  log transformed or normalized)}

  \item{exTable}{data frame, the first column contains
  names of genes or transcripts (Feature) and the remaining
  columns are counts for sample replicates spiked with ERCC
  controls}

  \item{repNormFactor}{optional vector of normalization
  factors for each replicate, default value is NULL and
  75th quantile normalization will be applied to
  replicates}

  \item{filenameRoot}{string root name for output files}

  \item{sample1Name}{string name for sample 1 in the gene
  expression experiment}

  \item{sample2Name}{string name for sample 2 in the gene
  expression experiment}

  \item{erccmix}{Name of ERCC mixture design, "RatioPair"
  is default, the other option is "Single"}

  \item{erccdilution}{unitless dilution factor used in
  dilution of the Ambion ERCC spike-in mixture solutions}

  \item{spikeVol}{volume in microliters of diluted ERCC mix
  spiked into the total RNA samples}

  \item{totalRNAmass}{mass in micrograms of total RNA
  spiked with diluted ERCC mixtures}

  \item{choseFDR}{False Discovery Rate for differential
  expression testing , default is 0.05}

  \item{ratioLim}{Limits for ratio axis on MA plot, default
  is c(-4,4)}

  \item{signalLim}{Limits for ratio axis on MA plot,
  default is c(-12,12)}

  \item{userMixFile}{optional filename input, default is
  NULL, if ERCC control ratio mixtures other than the
  Ambion product were used then a userMixFile can be used
  for the analysis}
}
\description{
Initialize the exDat list
}

