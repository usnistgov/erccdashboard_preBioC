dynRangePlotLODR <- function(dynRangeRes = dynRangeRes,LODR.annot.ERCC = LODR.annot.ERCC){
  folds = expDat$sampleInfo$FCcode
  legendLabels = expDat$sampleInfo$legendLabels
  lodrDat = LODR.annot.ERCC$LODRtable
  
#   myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
#   names(myColors) <- levels(FCcode$Ratio)
#   myColorsDiff <- myColors[-(which(FCcode$FC == 1))]
#   legendLabelsDiff <- legendLabels[-(which(FCcode$FC == 1))]
#   
#   colScale <- scale_colour_manual(name = "Ratio",values = myColorsDiff, labels = legendLabelsDiff)
#   fillScale <- scale_fill_manual(name = "Ratio", values = myColorsDiff, labels = legendLabelsDiff)
#   
  FCcodeC = folds[-c(which(folds$FC == 1)),]
  
  hDat = data.frame(Ratio = folds$Ratio, x = lodrDat$Log2Conc, y = lodrDat$Log2Count_normalized, xend = max(dynRangeRes$coordinates$limits$x), yend = lodrDat$Log2Count_normalized)
  #hDat = hDat[-which(hDat$FC == "1"),]
  hDat = hDat[which(is.finite(hDat$x)),]
  
  vDat = data.frame(Ratio = folds$Ratio, x = lodrDat$Log2Conc, y = min(dynRangeRes$coordinates$limits$y), xend = lodrDat$Log2Conc, yend = lodrDat$Log2Count_normalized )
  
  #vDat = vDat[-which(vDat$Ratio == "1"),]
  vDat = vDat[which(is.finite(vDat$x)),]
  
 if(dim(vDat)[1] == 0){ print("Error! Estimated distribution of p-values does not cross threshold p-value, may be due to insufficient data quantity"); break}
  if(dim(hDat)[1] == 0){ print("Error! Estimated distribution of p-values does not cross threshold p-value, may be due to insufficient data quantity"); break}
  
  dynRangeAnnot <- dynRangeRes + geom_segment(data = hDat, aes(x = x,y = y,xend = xend , yend = yend, colour = Ratio), lineend = "round",size = 2, alpha = 0.5) + geom_segment(data = vDat, aes(x = x,y = yend,xend = xend , yend = y, colour = Ratio),lineend = "round",arrow = arrow(length =unit(0.5,"cm")), size = 2,alpha = 0.5) + theme(legend.justification=c(0,1), legend.position=c(0,1))
  #print(dynRangeAnnot)
  expDat$Figures$plotdynRangeAnnot <- dynRangeAnnot
  return(expDat)
  
}
