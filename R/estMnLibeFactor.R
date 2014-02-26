estMnLibeFactor <- function( expDat, cnt = expDat$TranscriptsAB){
  
  sampleInfo = expDat$sampleInfo
  datType <- sampleInfo$datType
  
  if(datType == "count"){
    if(sampleInfo$totalSeqReads == F){
      #sampleLibeSums = colSums(cnt[-c(1)],na.rm =T)
      mnLibeFactor = (mean(as.vector(sampleLibeSums)))/10^6
      cat("\nUsing total mapped reads,\n mean library size factor = ")
      cat(mnLibeFactor)
    }else{
      mnLibeFactor = (mean(as.vector(expDat$totalReads)))/10^6
      cat("\nUsing total sequencing reads,\n mean library size factor = ")
      cat(mnLibeFactor)
      #sampleLibeSums = expDat$totalReads
    }  
  }
  if(datType == "array"){
    mnLibeFactor = mean(expDat$libeSize)
  }
  if(datType != "array"&&datType != "count")
    stop("Wrong datType, choose array or count")
  expDat$mnLibeFactor <- mnLibeFactor
  #expDat$sampleLibeSums <- sampleLibeSums
  return(expDat)
}
  
  