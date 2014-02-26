\name{runDashboard}
\alias{runDashboard}
\title{Run default erccdashboard analysis of ERCC control ratio mixtures}
\usage{
runDashboard(datType = NULL, expTable = NULL, repNormFactor = NULL,
  filenameRoot = NULL, sample1Name = NULL, sample2Name = NULL,
  erccmix = "RatioPair", erccdilution = 1, spikeVol = 1,
  totalRNAmass = 1, choseFDR = 0.05, userMixFile = NULL)
}
\arguments{
  \item{datType}{type is "count" or "array", unnormalized
  data is expected (normalized data may be accepted in
  future version of the package). Default is "count"
  (integer count data),"array" is unnormalized fluorescent
  intensities from microarray fluorescent intensities (not
  log transformed or normalized)}

  \item{expTable}{data frame, the first column contains
  names of genes or transcripts (Feature) and the remaining
  columns are counts for sample replicates spiked with ERCC
  controls}

  \item{repNormFactor}{vector of normalization factors for
  each replicate}

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
  expression testing}

  \item{userMixFile}{optional filename input, default is
  NULL, if ERCC control ratio mixtures other than the
  Ambion product were used then a userMixFile can be used
  for the analysis}
}
\description{
Run default erccdashboard analysis of ERCC control ratio
mixtures
}
\examples{
load(file = system.file("data/SEQC.Example.RData",
     package = "erccdashboard"))

expDat = runDashboard(datType = "count",
                 expTable = COH.RatTox.ILM.MET.CTL.countTable,
                 repNormFactor = COH.RatTox.ILM.MET.CTL.totalReads,
                 filenameRoot = "COH.ILM",
                 sample1Name = "MET", sample2Name = "CTL",
                 erccmix = "RatioPair", erccdilution = 1/100,
                 spikeVol = 1, totalRNAmass = 0.500,choseFDR = 0.1)

summary(expDat)
}

