erccdashboard
=============
This is the development version of the erccdashboard R package (v. 0.99.0). 
Use this package for analysis of the ERCC spike-in controls in differential 
gene expression experiments on any technology platform (including RNA-Seq or 
microarray experiments). The software automatically outputs performance 
measures derived from ERCC controls. These performance measures are output in 
a summary pdf file, results are also saved as Rdata for further analysis and to
enable "reproducible research"! These performance measures serve as an 
ERCC "dashboard" to enable scientists to understand technical performance of any 
differential gene expression experiment and to also compare experiments across
space and time.

The preprint of our manuscript describing this package is available on the arXiv
preprint server at this link: http://arxiv.org/abs/1406.4893

Please note that this is the development version of the code and it may change
frequently. Information will be posted here when the stable released version of
the package is available through an online R package repository.

Installation
------------
To install the erccdashboard locally in R:

1. Load devtools. You may need to run install.packages("devtools") 
if you don't already have devtools on your local machine:

    library("devtools")
    
2. Ensure that the three bioconductor packages that are erccdashboard
dependencies are installed using biocLite:
    
    source("http://bioconductor.org/biocLite.R")
    
    biocLite(c("edgeR","limma","qvalue"))

3. Use devtools to install the erccdashboard package development version:

    install_github("munrosa/erccdashboard")
    
All remaining erccdashboard package dependencies should be downloaded from CRAN.

Note that the installation may take some time, because of the examples that 
need to be run to build the user manual (vignette). Please be patient

4. Load the erccdashboard package
    
    library("erccdashboard")

5. Once the package is installed, the package vignette can be viewed, which 
provides detailed examples. Open the vignette pdf file with the command:
    
    vignette("erccdashboard")

Quick Start
----------
Reading the vignette is a good first step, but you can also use the following
steps as a quick start to test the package:

1. Load the example data

    data(SEQC.Example)

2. Run the dashboard demo on rat toxicogenomics data:

    exDat <- runDashboard(datType="count", isNorm=FALSE,
                           exTable=MET.CTL.countDat,
                           filenameRoot="RatTox", sample1Name="MET",
                           sample2Name="CTL", erccmix="RatioPair",
                           erccdilution=1/100, spikeVol=1,
                           totalRNAmass=0.500, choseFDR=0.1)

3. Review the results of the analysis in the RatTox.MET.CTL.RData and
   RatTox.MET.CTL.PDF file of the main diagnostic plots

4. Analysis of additional example data sets is described in the package vignette.