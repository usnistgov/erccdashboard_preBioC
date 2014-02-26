\name{erccdashboard-package}
\alias{erccdashboard-package}
\alias{erccdashboard}
\docType{package}
\title{
Assessment of Gene Expression Measurement Technical Performance with
External Spike-in RNA Controls
}
\description{
erccdashboard produces graphs and performance metrics for gene
    expression measurements derived from analysis of the External RNA
    Control Consortium (ERCC) spike-in controls.
}
\details{
\tabular{ll}{
Package: \tab erccdashboard\cr
Type: \tab Package\cr
Version: \tab 0.9.4\cr
URL: \tab http://nist.gov/erccControls\cr
Depends: \tab 
R (>= 3.0), ggplot2, gridExtra\cr
Imports: \tab 
reshape2,
plyr,
scales,
edgeR,
qvalue,
locfit,
QuasiSeq,
grid,
stringr,
ROCR,
gtools\cr
License: \tab GPL (>=2)\cr
Collate: \tab 
QL.results.R compareReps.R dashboardFile.R dynRangePlot.R
dynRangePlotLODR.R erccROC.R estLODR.R estMnLibeFactor.R est_r_m.R
geneExprTest.R getDesignMat.R initDat.R libeSizeNorm.R loadERCCInfo.R
loadExpMeas.R loadExpMeasFromRaw.R maConcPlot.R meltExpDat.R multiplot.R
prepERCCDat.R printLODRres.R saveInputs.R saveResults.R saveERCCPlots.R
sumTechReps.R testDECount.R testDEArray.R annotLODR.R plotAdjust.R
maSignalNew.R belowLODR.R runDashboard.R withinMixRatios.R
withinMixRatiosPairwise.R\cr
Built: \tab R 3.0.2; ; 2014-02-14 21:36:51 UTC; unix\cr
}

Index:
\preformatted{
ERCCDef                 ERCC lengths and GC content
ERCCMix1and2            Ambion ERCC control ratio mixture design and
                        transcript concentrations
annotLODR               Annotate signal-abundance and ratio-abundance
                        plots with LODR
dynRangePlot            Produce signal-abundance plot to evaluate
                        dynamic range
erccROC                 Produce Receiver Operator Characteristic (ROC)
                        Curves and AUC statistics
estLODR                 Estimate Limit of Detection of Ratios (LODR)
est_r_m                 Estimate the mRNA fraction differences for the
                        pair of samples using replicate data
geneExprTest            Prepare differential expression testing results
                        for spike-in analysis
initDat                 Initialize the expDat list
maSignal                Generate MA plots with or without annotation
                        using LODR estimates
runDashboard            Run default erccdashboard analysis of ERCC
                        control ratio mixtures
saveERCCPlots           Save erccdashboard plots to a pdf file
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{erccdashboardVignette} \tab erccdashboard examples (source, pdf)\cr
}
}
\author{
Sarah Munro, Steve Lund

Maintainer: Sarah Munro <sarah.munro@nist.gov>
}
\references{
https://www-s.nist.gov/srmors/certificates/2374.pdf
}