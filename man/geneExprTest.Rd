% Generated by roxygen2 (4.0.2): do not edit by hand
\name{geneExprTest}
\alias{geneExprTest}
\title{Prepare differential expression testing results for spike-in analysis}
\usage{
geneExprTest(exDat)
}
\arguments{
\item{exDat}{list, contains input data and stores analysis results}
}
\description{
Prepare differential expression testing results for spike-in analysis
}
\details{
This function wraps the QuasiSeq differential expression testing package for
datType = "count" or uses the limma package for differential expression
testing if datType = "array". Alternatively, for count data only, if
correctly formatted DE test results are provided,
then geneExprTest will bypass DE testing (with reduced runtime).
}
\examples{
data(SEQC.Example)

exDat <- initDat(datType="array", isNorm=FALSE,
                 exTable=UHRR.HBRR.arrayDat,
                 filenameRoot="testRun", sample1Name="UHRR",
                 sample2Name="HBRR", erccmix="RatioPair",
                 erccdilution = 1, spikeVol = 50,
                 totalRNAmass = 2.5*10^(3), choseFDR=0.01)

exDat <- est_r_m(exDat)

exDat <- dynRangePlot(exDat)

exDat <- geneExprTest(exDat)
}

