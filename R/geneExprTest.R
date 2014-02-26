#' Prepare differential expression testing results for spike-in analysis
#'
#' @param expDat    list, contains input data and stores analysis results
#' 
#' @details
#' This function wraps the QuasiSeq differential expression testing package for
#' datType = "count" or uses the limma package for differential expression 
#' testing if datType = "array". Alternatively, for count data only, if a 
#' correctly formatted csv file is provided with the necessary DE test results, 
#' then geneExprTest will bypass DE testing (with reduced runtime). The function
#' will look for a csv file with the name "filenameRoot.quasiSeq.res.csv" and 
#' columns corresponding to "Feature", "pvals", and "qvals" should be in the 
#' file.
#' 
#' @export

geneExprTest <- function(expDat){
  datType <- expDat$sampleInfo$datType
  if(datType == "array"){
    if(is.null(expDat$sampleInfo$choseFDR)){
      getPThresh<- function(){
        cat("\nFDR is NULL, to continue with P-values for LODR estimation\n")
        readline("Enter the threshold P-value: ")
      }
      expDat$Results$p.thresh <- as.numeric(getPThresh())  
    }
    
    expDat <- testDEArray(expDat)
    cat("\nFinished DE testing\n")
  }
  if(datType == "count"){
    cnt <- expDat$Transcripts
    designMat <- expDat$designMat
    sampleInfo <- expDat$sampleInfo
    info <- designMat
    choseFDR <- sampleInfo$choseFDR
    
    # set initial p.thresh
    p.thresh<-.1
    # First 3 columns of qvalFile must contain Feature, pvals, and qvals
    
    qvalFile <- paste(sampleInfo$filenameRoot,"quasiSeq.res.csv",sep=".")
    pvalERCC <- paste(sampleInfo$filenameRoot, "ERCC Pvals.csv",sep=" ")
    #dispPlotFile = paste(sampleInfo$filenameRoot, ".DispPlot.RData",sep = "")
    
    # Decide to reuse results or run testDE
    if (file.exists(qvalFile) == TRUE){
      deRes <- read.csv(qvalFile)
      stopifnot("qvals" %in% names(deRes))
      if (file.exists(pvalERCC) == TRUE){
        cat(paste("\n Differential expression test results exist, will use  \n",
                  "existing P-values and Q-values for analysis.\n",
                  "Delete", qvalFile, "if you want to repeat differential \n",
                  "differential expression testing or view dispersion plots\n"))
        if(any(deRes$qvals<choseFDR)){
          p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
        }
      }else{
        cat("\nStarting differential expression tests\n")
        expDat <- suppressWarnings(testDECount(sampleInfo, expDat, cnt = cnt, 
                                          info = info ))  
        deRes <- read.csv(qvalFile)
        if(any(deRes$qvals<choseFDR)){
          p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
        }
      }
      
    }
    if(file.exists(qvalFile) == FALSE){
      pvalERCC = paste(sampleInfo$filenameRoot, "ERCC Pvals.csv",sep=" ")
      if(file.exists(pvalERCC)){
        cat(paste("\n P-values exist for ERCCs, but Q-values are missing\n"))
        # First three columns must be Feature, MnCnt, Pval
        # Need to get Threshold p-value from user
        getPThresh<- function(){
          cat("\nTo continue with P-values for LODR estimation\n")
          readline("Enter the threshold P-value: ")
        }
        p.thresh <- as.numeric(getPThresh())
        
      }else{
        cat("\nStarting differential expression tests\n")
        expDat <- suppressWarnings(testDECount(sampleInfo, expDat, cnt = cnt, 
                                          info = info ))  
        deRes <- read.csv(qvalFile)
        if(any(deRes$qvals<choseFDR)){
          p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
        } 
      }
      
    }
    
    
    if(is.null(expDat$Figures$dispPlot)){
      cat(paste("\nDE testing results supplied without companion dispersion\n",
                "plot. Dispersion plot is unavailable to print.\n"))
    }  
    cat("\nThreshold P-value\n")
    cat(p.thresh,"\n")
    
    if (p.thresh > .1){
      cat(paste("Threshold P-value is high for the chosen FDR of ", 
                as.character(choseFDR)))
      cat(paste("\nThe sample comparison indicates a large amount of \n",
                "differential expression in the measured transcript \n",
                "populations\n"))
    }
    expDat$Results$p.thresh <- p.thresh
  }
  
  return(expDat)
  
}
