#' Annotate signal-abundance and ratio-abundance plots with LODR
#'
#' @param expDat    list, contains input data and stores analysis results
#' 
#' @export
#' 
annotLODR <- function(expDat){
  filenamePval <- paste0(expDat$sampleInfo$filenameRoot,".quasiSeq.res.csv")
  cat("\n ERCC LODR estimates are available\n")
  LODR.annot.ERCC <- printLODRres(expDat)
 # expDat <- dynRangePlotLODR(dynRangeRes = expDat$Figures$plotdynRange,
 #                   LODR.annot.ERCC = LODR.annot.ERCC)
  
 # expDat<- maConcPlot(expDat, LODR.annot.ERCC, alphaPoint = 0.8, r_mAdjust = T, 
 #                       replicate = T)
  expDat$Results$LODR.annot.ERCC <- LODR.annot.ERCC
  
  expDat <- maSignal(expDat)
  
  expDat <- belowLODR(expDat,pvalDat=filenamePval)
  
  #print(expDat$Figures$maPlot)
  return(expDat)
}