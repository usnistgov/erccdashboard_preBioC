# erccdashboard

## About

Use the erccdashboard for analysis of ERCC spike-in controls in differential 
gene expression experiments on any technology platform (including RNA-Seq or 
microarray experiments). The software automatically outputs performance 
measures derived from ERCC controls to a summary pdf file and the 
results are also saved in an Rdata file for further analysis (and to enable 
"reproducible research"). These ERCC-derived performance measures serve as an 
ERCC "dashboard" to enable scientists to understand technical performance of any 
differential gene expression experiment and to also compare experiments across
space and time. When you use this package please cite our publication:

Munro, S. A. et al. Assessing technical performance in differential gene 
expression experiments with external spike-in RNA control ratio mixtures. 
Nat. Commun. 5:5125 doi: 10.1038/ncomms6125 (2014).

## Installation


### To install this GitHub development version: 
   
1. Get Bioconductor package dependencies
    ```    
    source("http://bioconductor.org/biocLite.R")
    
    biocLite(c("edgeR","limma","qvalue"))
    ```
2. Use the devtools command to install erccdashboard from github
    ```
    install_github("usnistgov/erccdashboard")
    ```

### To install the release version from Bioconductor:

1. Run these Bioconductor installation commands (all Bioconductor and CRAN 
dependencies will be automatically installed)
    ```
    source("http://bioconductor.org/biocLite.R")
    
    biocLite("erccdashboard")
    ```

## Quick Start

1. Load the erccdashboard package
    ```
    library("erccdashboard")
    ```
2. Load the example data
    ```
    data(SEQC.Example)
    ```
3. Run the dashboard demo on rat toxicogenomics data:
    ```
    exDat <- runDashboard(datType="count", isNorm=FALSE,
                           exTable=MET.CTL.countDat,
                           filenameRoot="RatTox", sample1Name="MET",
                           sample2Name="CTL", erccmix="RatioPair",
                           erccdilution=1/100, spikeVol=1,
                           totalRNAmass=0.500, choseFDR=0.1)
    ```
4. Review the results of the analysis in the RatTox.MET.CTL.RData and
   RatTox.MET.CTL.pdf file of the main diagnostic plots

5. Analysis of additional example data sets is described in the package vignette.
Open the vignette pdf file with the command:
    ```
    vignette("erccdashboard")
    ```