erccdashboard
=============
This is the source code for the erccdashboard R packge (v. 0.9.7)

Note that the erccdashboard package has other packages as dependencies. These packages will need to be installed and loaded by the user. Once the erccdashboard package is on a remote repository (e.g. CRAN or Bioconductor) the package dependencies will automatically be addressed during installation and loading of the erccdashboard package. All package dependencies are on CRAN except for edgeR and qvalue which are both on Bioconductor. Use the following R commands to install edgeR and qvalue: 

source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("qvalue")

The full list of packages that must be installed and loaded are:
ggplot2
reshape2
plyr
scales
edgeR
locfit
QuasiSeq
grid
gridExtra
stringr
ROCR
gtools
qvalue

Packages can be installed and loaded with the use of install.package() and library() functions.

Instructions for erccdashboard Installation and Quick Start:

Installation:

1. Download the source zip archive to your system.

2. At R prompt install the package with:

  > install.packages(pkgs="path_to_package/erccdashboard_0.9.7")

3. Load the package with

	> library(erccdashboard)

   Errors will occur if any of the above package dependencies are missing.

4. You can review the package vignette for an explanation of how to run the dashboard or follow this quick start:

Quick Start:

1. Load the example data

	> data(SEQC.Example)

2. Run the dashboard demo on rat toxicogenomics data:

	> expDat <- runDashboard(datType = "count",
                       expTable = COH.RatTox.ILM.MET.CTL.countTable,
                       filenameRoot = "COH.ILM",sample1Name = "MET",
                       sample2Name = "CTL",erccmix = "RatioPair",
                       erccdilution = 1/100, spikeVol = 1,
                       totalRNAmass = 0.500, choseFDR = 0.1)

3. Review the results of the analysis in the COH.ILM.MET.CTL.RData and
   COH.ILM.MET.CTL.PDF file of the main diagnostic plots

4. Analysis of different example data sets is described in the package vignette.