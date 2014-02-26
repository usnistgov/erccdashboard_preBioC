#' Save erccdashboard plots to a pdf file
#'
#' @param expDat     list, contains input data and stores analysis results
#' @param plotsPerPg string, if "manuscript" then the 4 main plots are printed  
#'                   to one page, if "single" then each plot is printed to page 
#'                   in the pdf file                     
#' @param plotlist   list, contains plots to print
#' 
#' @description
#' The function savePlots will save selected figures to a pdf file. The default 
#' is the 4 manuscript figures to a single page (plotsPerPg = "manuscript"). 
#' If plotsPerPg = "single" then each plot is placed on an 
#' individual page in one pdf file. If plotlist is not defined (plotlist = NULL)
#'  then all plots in expDat$Figures are printed to the file.
#' 
#' # to print 4 plots from manuscript to a single page pdf file
#' saveERCCPlots(expDat, plotsPerPg = "manuscript")
#' 
#' # to create a multiple page pdf of all plots produced
#' saveERCCPlots(expDat, plotsPerPg = "single", plotlist = expDat$Figures)
#' 
#' @export

saveERCCPlots<-function(expDat,plotsPerPg = "manuscript", plotlist = NULL){
  #Options are either the default of printing the plots as shown in publication
  # plotsPerPg = "manuscript" and plotlist is NULL or plotsPerPg = "single" and
  # any combination of the plots can be printed, one per page
  # Open PDF file to write results
  filenameUse <- expDat$sampleInfo$filenameRoot 
#   if (plotsPerPg == "manuscript"){
#     cols = 2
#     pwidth = 7*cols
#     pheight = 7*6/cols
#     pdf(file = paste(filenameUse,"pdf",sep="."),title=filenameUse, 
#         width=pwidth,height = pheight)
#     
#     multiplot(expDat$Figures$rocPlot,expDat$Figures$dynRangePlot, 
#               expDat$Figures$lodrERCCPlot,expDat$Figures$rangeResidPlot, 
#               expDat$Figures$dispPlot,expDat$Figures$maPlot,cols=2)
#     dev.off()
#   } 

cat("\nSaving main dashboard plots to pdf file...")
  if (plotsPerPg == "manuscript"){
    cols = 2
    nFigs = 4
    pwidth = 7*cols
    pheight = 7*nFigs/cols
    pdf(file = paste(filenameUse,"pdf",sep="."),title=filenameUse,
        width=pwidth,height = pheight)
    #pdf(file =  paste(filenameUse,"pdf",sep="."),title=filenameUse,
    #    paper = "letter")
    multiplot(expDat$Figures$dynRangePlot, expDat$Figures$rocPlot,
              expDat$Figures$maPlot, expDat$Figures$lodrERCCPlot, cols=cols)
    dev.off()
  }
  if (plotsPerPg == "single"){
    if (is.null(plotlist)){
      plotlist = expDat$Figures
    } 
    pdf(file = paste(filenameUse,"pdf",sep="."),onefile=T,width=7,height = 7)
    print(plotlist)
    dev.off()
  }
  
}
