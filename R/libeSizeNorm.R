libeSizeNorm <- function(expDat){
  
  sampleInfo <- expDat$sampleInfo
  erccInfo <- expDat$erccInfo
  expressDat = expDat$Transcripts
  datType <- sampleInfo$datType
  repNormFactor <- sampleInfo$repNormFactor
  #print(sampleInfo$repNormFactor)
  normVec = T
  if(sampleInfo$datType == "count"){
    if (is.null(sampleInfo$repNormFactor)){
      normVec = F
    }
   # print(normVec)
  }
  if(sampleInfo$datType == "array"){
    normVec = F
  }
  
  # Library size normalize the data
  if (sampleInfo$libeSizeNorm == T){
    if(sampleInfo$datType == "array"){
      #cat("\nUsing median intensity for ERCC 1:1 controls\n")
      #cat("\nUse median intensity from each array for normalization\n")
      #cat("\nUsing total intensity to normalize each array\n")
      cat("\nUsing 75th percentile (upper quartile) intensity to normalize each array\n")
      ### Subset the 1:1 ercc controls and calculate the median intensity for
      # for each column, build a libeSize vector  
      TranscriptsAll = expressDat
      #ercc1to1 <- expressDat[which(erccInfo$idColsSRM$Ratio == "1:1"),]
      #libeSize = apply(ercc1to1[-c(1)],MARGIN=2,FUN=median)
      libeSize = apply(expressDat[-c(1)],MARGIN=2,FUN=quantile,probs=0.75)
      #libeSize  = colSums(TranscriptsAll[-c(1)])
      #print(libeSize)
      
      libAdjust = sweep(expressDat[-c(1)],2,libeSize,"/")
      expressDat = cbind(expressDat[c(1)], libAdjust)
    }
    if(sampleInfo$datType == "count"){
      if (normVec == F){
        cat(paste("\nrepNormFactor is NULL,\n",
                  "Using Default Upper Quartile Normalization Method",
                  " - 75th percentile)\n"))
        TranscriptsAll = expressDat
        libeSize = apply(expressDat[-c(1)],MARGIN=2,FUN=quantile,probs=0.75)
        #TranscriptMappedReadSums = colSums(TranscriptsAll[-c(1)],na.rm = T)
        #libeSize = TranscriptMappedReadSums
        datCols = expressDat[-c(1)]
        libeSize = libeSize#/(10^6) #per million mapped reads
        #Library size normalize the data  
        libAdjust = sweep(datCols, 2, libeSize,"/")
        expressDat = cbind(expressDat[c(1)], libAdjust)
      }
      
      if (normVec == T){
        cat("\nUsing read depth normalization factors provided in repNormFactor\n")
        TranscriptsAll = expressDat
        libeSize = repNormFactor
        #print(libeSize)
        datCols = expressDat[-c(1)]
        libeSize = libeSize/(10^6) #per million mapped reads
        
        #Library size normalize the data  
        libAdjust = sweep(datCols, 2, libeSize,"/")
        expressDat = cbind(expressDat[c(1)], libAdjust)
      }
    }
    
    cat("\nLibrary sizes:\n")
    cat(libeSize)
    
  }
  else{
    libeSize = NULL
  }
    # get just ERCC data in the expression data frame
  expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]
  
  
  mnLibeFactor = mean(libeSize)
  
  expDat$normERCCDat <- expressDat
  expDat$libeSize <- libeSize
  expDat$mnLibeFactor <- mnLibeFactor

  
  return(expDat)
  
}