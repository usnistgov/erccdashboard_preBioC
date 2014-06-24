erccdashboard
=============
This is the development version of the erccdashboard R package (v. 0.9.10).
Information will be posted here when the stable released version of the package
is available.
Feel free to use the development version of the code.

Installation
------------
To install the erccdashboard locally in R:

1. Load devtools. You may need to run install.packages("devtools"") 
if you don't already have devtools on your local machine:

    library("devtools")

2. Use devtools to install this development version:

    install_github("munrosa/erccdashboard")


Note that the installation may take some time, because of the examples that 
need to be run to build the user manual (vignette). Please be patient.

3. Once the package is installed, the package vignette can be viewed, which 
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
                           exTable=MET.CTL.countTable,
                           filenameRoot="RatTox", sample1Name="MET",
                           sample2Name="CTL", erccmix="RatioPair",
                           erccdilution=1/100, spikeVol=1,
                           totalRNAmass=0.500, choseFDR=0.1)

3. Review the results of the analysis in the RatTox.MET.CTL.RData and
   RatTox.MET.CTL.PDF file of the main diagnostic plots

4. Analysis of additional example data sets is described in the package vignette.