dashboardFile <-function(expDat, filenameRoot){
  sampleInfo <- expDat$sampleInfo
 
  filenameUse = with(sampleInfo, paste(filenameRoot,sample1Name,
                                          sample2Name,sep = "."))  
  cat(paste("Filename root is:", filenameUse, "\n"))
  
  expDat$sampleInfo$filenameRoot <- filenameUse
  
  return(expDat)
}