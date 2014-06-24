#' Annotate signal-abundance and ratio-abundance plots with LODR
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @export
#' 
annotLODR <- function(exDat){
    
    ## Assign LODR results to object in exDat
    LODR.annot.ERCC <- printLODRres(exDat)
    
    exDat$Results$LODR.annot.ERCC <- LODR.annot.ERCC
    
    ## Produce MA plots
    exDat <- maSignal(exDat)
    
    return(exDat)
}