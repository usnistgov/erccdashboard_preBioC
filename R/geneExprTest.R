#' Prepare differential expression testing results for spike-in analysis
#'
#' @param exDat    list, contains input data and stores analysis results
#' 
#' @details
#' This function wraps the QuasiSeq differential expression testing package for
#' datType = "count" or uses the limma package for differential expression 
#' testing if datType = "array". Alternatively, for count data only, if
#' correctly formatted DE test results are provided,
#' then geneExprTest will bypass DE testing (with reduced runtime).
#' 
#' @export

geneExprTest <- function(exDat){
  datType <- exDat$sampleInfo$datType
  isNorm <- exDat$sampleInfo$isNorm
  #print(datType)
  if(datType == "array"){
    if(is.null(exDat$sampleInfo$choseFDR)){
      getPThresh<- function(){
        cat("\nFDR is NULL, to continue with P-values for LODR estimation\n")
        readline("Enter the threshold P-value: ")
      }
      exDat$Results$p.thresh <- as.numeric(getPThresh())  
    }
    
    exDat <- testDEArray(exDat)
    cat("\nFinished DE testing\n")
  }
  if(datType == "count"){
    cnt <- exDat$Transcripts
    designMat <- exDat$designMat
    sampleInfo <- exDat$sampleInfo
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
                  "expression testing or view dispersion plots\n"))
        if(any(deRes$qvals<choseFDR)){
          p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
        }
      }else{
        cat("\nStarting differential expression tests\n")
        exDat <- suppressWarnings(testDECount(sampleInfo, exDat, cnt = cnt, 
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
        cat(paste("\n Using existing ERCC DE test p-values for AUC and LODR\n"))
        # First three columns must be Feature, MnSignal, Pval
        pvalFile <- read.csv(pvalERCC)
        pvalFile$qvals <- qvalue(pvalFile$Pval)$qvalues
        if(any(pvalFile$qvals<choseFDR)){
          p.thresh<-max(pvalFile$Pval[pvalFile$qvals<choseFDR])
        }
        
        # Need to get Threshold p-value from user
        #getPThresh<- function(){
        #  cat("\nTo continue with P-values for LODR estimation\n")
        #  readline("Enter the threshold P-value: ")
        #}
        #p.thresh <- as.numeric(getPThresh())
        
      }else{
        if (isNorm == T){
          cat(paste0("\nQuasiSeq DE Testing requires count (integer) data.\n",
                     "To estimate AUC and LODR for normalized RNA-Seq data\n",
                     "the file '",sampleInfo$filenameRoot,
                     " ERCC Pvals.csv' is required\n"))
          return(exDat)
        }
        cat("\nStarting differential expression tests\n")
        exDat <- suppressWarnings(testDECount(sampleInfo, exDat, cnt = cnt, 
                                          info = info ))  
        deRes <- read.csv(qvalFile)
        if(any(deRes$qvals<choseFDR)){
          p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
        } 
      }
      
    }
    
    
    if(is.null(exDat$Figures$dispPlot)){
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
    exDat$Results$p.thresh <- p.thresh
  }
  
  return(exDat)
  
}
