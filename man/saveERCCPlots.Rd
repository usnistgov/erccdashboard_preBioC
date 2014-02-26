\name{saveERCCPlots}
\alias{saveERCCPlots}
\title{Save erccdashboard plots to a pdf file}
\usage{
saveERCCPlots(expDat, plotsPerPg = "manuscript", plotlist = NULL)
}
\arguments{
  \item{expDat}{list, contains input data and stores
  analysis results}

  \item{plotsPerPg}{string, if "manuscript" then the 4 main
  plots are printed to one page, if "single" then each plot
  is printed to page in the pdf file}

  \item{plotlist}{list, contains plots to print}
}
\description{
The function savePlots will save selected figures to a pdf
file. The default is the 4 manuscript figures to a single
page (plotsPerPg = "manuscript"). If plotsPerPg = "single"
then each plot is placed on an individual page in one pdf
file. If plotlist is not defined (plotlist = NULL) then all
plots in expDat$Figures are printed to the file.

# to print 4 plots from manuscript to a single page pdf
file saveERCCPlots(expDat, plotsPerPg = "manuscript")

# to create a multiple page pdf of all plots produced
saveERCCPlots(expDat, plotsPerPg = "single", plotlist =
expDat$Figures)
}

