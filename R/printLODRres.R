printLODRres <- function(exDat){
  sampleInfo <- exDat$sampleInfo
  erccInfo <- exDat$erccInfo
  
  fit.coeff <- exDat$fit.coeff
  mnLibeFactor <- exDat$mnLibeFactor
  
  FCcode = erccInfo$FCcode
  legendLabels = erccInfo$legendLabels
  

  lodr.res = data.frame(exDat$Results$lodr.res.ERCC)
  
  ### Fold = lodr.res[c(1)]
  Fold = as.numeric(exDat$erccInfo$FCcode$FC)

  Count = as.numeric(gsub("<", "",lodr.res$Estimate))
  Ratio = legendLabels 
  
  
  # Convert LODR count estimate to library size normalized
  logCount = log2((Count/(mnLibeFactor)))#+1)
  
  ###ConcEst = (log2(Count/(mnLibeFactor))-fit.coeff[1])/fit.coeff[2]
  
  ###LODR.print.res = data.frame(Fold, Ratio, Count, logCount, ConcEst)
  LODR.print.res = data.frame(Fold, Ratio, Count, logCount)
  
  #names(LODR.print.res)<- c("Fold","Ratio","Count","Log2Count_normalized","Log2Conc")
  names(LODR.print.res)<- c("Fold","Ratio","Count","Log2Count_normalized")
  print(LODR.print.res)
  
  #cutoffs = ConcEst[which(!(is.na(ConcEst)))]
  countCutoffs = logCount[which(!(is.na(logCount)))]
  return(list(LODRtable = LODR.print.res, 
              #cutoffs=cutoffs,
              countCutoffs = countCutoffs))
  
}