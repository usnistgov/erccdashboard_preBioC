#' Annotate signal-abundance and ratio-abundance plots with LODR
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @export
#' 
annotLODR <- function(exDat){
  #filenamePval <- paste0(exDat$sampleInfo$filenameRoot,".quasiSeq.res.csv")
  
  LODR.annot.ERCC <- printLODRres(exDat)

  exDat$Results$LODR.annot.ERCC <- LODR.annot.ERCC

  exDat <- maSignal(exDat)
  
  return(exDat)
}