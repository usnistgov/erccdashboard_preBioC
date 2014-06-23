dashboardFile <-function(expDat, filenameRoot){
  sampleInfo <- expDat$sampleInfo

  filenameUse = paste(filenameRoot,sampleInfo$sample1Name,
                      sampleInfo$sample2Name,sep = ".")  
  cat(paste("Filename root is:", filenameUse, "\n"))
  
  expDat$sampleInfo$filenameRoot <- filenameUse
  
  return(expDat)
}