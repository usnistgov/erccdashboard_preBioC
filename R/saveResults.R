# Save erccdashboard analysis results to a RData file for future analysis
saveResults <- function(expDat){
  
  filenameRoot <- expDat$sampleInfo$filenameRoot
  # Name and consolidate metrics for the interlaboratory comparison
  #nam <- paste(filenameRoot, "expDat",sep = ".")
  #assign(nam,expDat)
  
  nam <- paste(filenameRoot, "AUC",sep = ".")
  assign(nam,expDat$AUCdat)
  
  nam <- paste(filenameRoot, "lodr.ERCC",sep = ".")
  assign(nam,expDat$lodr.res.ERCC)
  
  nam <- paste(filenameRoot, "sdGlobal",sep = ".")
  assign(nam,expDat$sdGlobal)
  
  nam <- paste(filenameRoot, "modRatVar",sep = ".")
  assign(nam,expDat$modRatVar)
  
  nam <- paste(filenameRoot, "r_m",sep = ".")
  assign(nam,expDat$r_m.res$r_m.mn)
  
  nam <- paste(filenameRoot, "r_m_lower",sep = ".")
  assign(nam,expDat$r_m.res$r_m.lower)
  
  nam <- paste(filenameRoot, "r_m_upper",sep = ".")
  assign(nam,expDat$r_m.res$r_m.upper)
  
  nam <- paste(filenameRoot, "p.thresh",sep = ".")
  assign(nam,expDat$p.thresh)
  
  nam <- paste(filenameRoot, "figures",sep = ".")
  assign(nam,expDat$Figures)
  
  to.save <- ls()
  
  save(list = to.save[grepl(pattern = filenameRoot,x=to.save)],
       file=paste(filenameRoot,"Results","RData",sep = "."))
  
}