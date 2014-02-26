saveInputs <- function(expDat){
  study <- expDat$sampleInfo$study
  filenameRoot <- expDat$sampleInfo$filenameRoot
  
  # Name and consolidate metrics for the interlaboratory comparison
  #nam <- paste(filenameRoot, "expDat",sep = ".")
  #assign(nam,expDat)
  
  if (study == "SEQC_Main"){
    if (dim(expDat$expressDatSumNoNorm[-c(1)])[2] != length(expDat$totalReadSum)){
      print("countTable and totalReads array do not match!")
      break
    }
    
    nam <- paste(filenameRoot, "countTable",sep = ".")
    assign(nam,expDat$expressDatSumNoNorm)
    
    nam <- paste(filenameRoot, "totalReads",sep = ".")
    assign(nam,expDat$totalReadSum)
    
  }
  if (study == "SEQC_RatTox"){
    # To Do Sarah Munro edit with correct parts of expDat for RatTox data
    if (dim(expDat$TranscriptsAB[-c(1)])[2] != length(expDat$totalReads)){
      print("countTable and totalReads array do not match!")
      break
    }
    
    nam <- paste(filenameRoot, "countTable",sep = ".")
    assign(nam,expDat$TranscriptsAB)
    
    nam <- paste(filenameRoot, "totalReads",sep = ".")
    assign(nam,expDat$totalReads)
    
  }
  
  to.save<- ls()
  save(list = to.save[grepl(pattern = filenameRoot,x=to.save)],
       file=paste(filenameRoot,"Inputs","RData",sep = "."))
    
}