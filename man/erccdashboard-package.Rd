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
Version: \tab 0.9.7\cr
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
compareReps.R dashboardFile.R dynRangePlot.R
dynRangePlotLODR.R erccROC.R estLODR.R estMnLibeFactor.R est_r_m.R
geneExprTest.R getDesignMat.R initDat.R libeSizeNorm.R loadERCCInfo.R
loadExpMeas.R loadExpMeasFromRaw.R maConcPlot.R meltExpDat.R multiplot.R
prepERCCDat.R printLODRres.R saveInputs.R saveResults.R saveERCCPlots.R
sumTechReps.R testDECount.R testDEArray.R annotLODR.R plotAdjust.R
maSignalNew.R belowLODR.R runDashboard.R withinMixRatios.R
withinMixRatiosPairwise.R quasiUpdateQLResults.R\cr
Built: \tab R 3.0.2; ; 2014-04-25 21:24:33 UTC; unix\cr
}

Index:
\preformatted{
annotLODR               Annotate signal-abundance and ratio-abundance
                        plots with LODR
dynRangePlot            Produce signal-abundance plot to evaluate
                        dynamic range
erccROC                 Produce Receiver Operator Characteristic (ROC)
                        Curves and AUC statistics
erccdashboard-package   Assessment of Gene Expression Measurement
                        Technical Performance with External Spike-in
                        RNA Controls
estLODR                 Estimate Limit of Detection of Ratios (LODR)
est_r_m                 Estimate the mRNA fraction differences for the
                        pair of samples using replicate data
geneExprTest            Prepare differential expression testing results
                        for spike-in analysis
initDat                 Initialize the exDat list
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

To cite the erccdashboard analysis package:\cr
Munro SA, Lund S, Pine PS, Binder H, Clevert D, Conesa A, Dopazo J, Fasold
  M, Hochreiter S, Hong H, Jafari N, Kreil DP, Łabaj PP, Liao Y, Lin S,
  Meehan J, Mason CE, Santoyo J, Setterquist RA, Shi L, Shi W, Smyth GK,
  Stralis-Pavese N, Su Z, Tong W, Wang C, Wang J, Xu J, Ye Z, Yang Y, Yu Y,
  Salit M (Under Review, 2014). Assessing Technical Performance in Gene
  Expression Experiments with External Spike-in RNA Control Ratio Mixtures.\cr

The certificate for the NIST Standard Reference Material (SRM 2374) which is the
plasmid DNA library source of the ERCC control sequences is available here:\cr
https://www-s.nist.gov/srmors/certificates/2374.pdf\cr
  
}