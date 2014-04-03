#' Annotate signal-abundance and ratio-abundance plots with LODR
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @export
#' 
annotLODR <- function(exDat){
  filenamePval <- paste0(exDat$sampleInfo$filenameRoot,".quasiSeq.res.csv")
  cat("\n ERCC LODR estimates are available\n")
  LODR.annot.ERCC <- printLODRres(exDat)
 # exDat <- dynRangePlotLODR(dynRangeRes = exDat$Figures$plotdynRange,
 #                   LODR.annot.ERCC = LODR.annot.ERCC)
  
 # exDat<- maConcPlot(exDat, LODR.annot.ERCC, alphaPoint = 0.8, r_mAdjust = T, 
 #                       replicate = T)
  exDat$Results$LODR.annot.ERCC <- LODR.annot.ERCC
  
  exDat <- maSignal(exDat)
  
  #exDat <- belowLODR(exDat,pvalDat=filenamePval)
  
  #print(exDat$Figures$maPlot)
  return(exDat)
}