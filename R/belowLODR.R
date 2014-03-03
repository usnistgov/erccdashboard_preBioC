belowLODR <- function(expDat,pvalDat = NULL){
  
  datType<- expDat$sampleInfo$datType
  p.thresh <- expDat$Results$p.thresh
  choseFDR <- expDat$sampleInfo$choseFDR
  LODRtable <- expDat$Results$LODR.annot.ERCC$LODRtable
  # Remove Fold change of 1 (Ratio is 1:1)
  LODRtable<- LODRtable[(LODRtable$Fold!=1),]

  Ratio <- as.character(LODRtable$Ratio)
  FC <-LODRtable$Fold
  
  if (datType == "count"){
    if(missing(pvalDat))stop("Provide filename as pvalDat argument!")
    pval <- read.csv(pvalDat, header = T)
    pval$Log2Rat <- abs(pval$Log2Rat)
  }
  if (datType == "array"){
    pval = expDat$Results$limma.res
    pval$A <- 2^(expDat$Results$limma.res$AveExpr)
    pval$logFC <- abs(pval$logFC)
  }
  #cat(paste("\nTotal Number of features tested for differential expression",
  #          "=", dim(pval)[1]))
  pvalCutoff <- NULL
  if(datType == "array"){
    pvalCutoff <- pval[(pval$adj.P.Val <= p.thresh),]
  }
  if(datType == "count"){
    pvalCutoff <- pval[(pval$pvals <= p.thresh),] 
  }
  pvalCutoff$Feature <- as.character(pvalCutoff$Feature)
  
  cat(paste("\nNumber of DE genes =",dim(pvalCutoff)[1] ,"(based on p.thresh =", 
            round(p.thresh,digits = 4), "for FDR =",choseFDR,") out of\n", dim(pval)[1],
            "total features tested for differential expression\n"))

  cat(paste("\nBegin flagging DE features using fold change, p.thresh,",
            "and LODR...\n"))
  cat("\n")
  for (i in 1:length(FC)){
    cat(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
    logFCcut <- abs(log2(FC[i]))
    #print(logFC)
    #estIdx <- match(logFCcut, abs(log2(expDat$erccInfo$FCcode$FC)))
    estLODRraw <- as.numeric(as.character(LODRtable$Count[i]))
    #estLODRnorm <- as.numeric(LODRtable$Log2Count_normalized[i])
    if (estLODRraw != "Inf"){
      cat(paste0("\nFor abs(log2(fold change)) less than or equal to ", 
                 round(logFCcut,digits=3), 
                " (Ratio = ",Ratio[i],"), p.thresh = ",
                round(p.thresh,digits=4),",\n"))
      
      
      if(datType == "array"){
        
        cat(paste("and below LODR estimate =",
                   estLODRraw,"average log2 fluorescent signal\n"))
        
        belowLODR <- pvalCutoff[((pvalCutoff$A < estLODRraw)&
                                   ((abs(pvalCutoff$logFC) < logFCcut)|
                                      (is.na(abs(pvalCutoff$logFC))))),]
#                             (logFC >= logFCcut) & 
#                               (A <= estLODRnorm ) & 
#                               (adj.P.Val <= p.thresh))
      }
      if(datType == "count"){
        
        cat(paste("and below LODR estimate =",
                  estLODRraw,"average counts"))
        
        belowLODR <- pvalCutoff[((pvalCutoff$MnCnt < estLODRraw)&
                                   ((abs(pvalCutoff$Log2Rat) < logFCcut)|
                                      (is.na(abs(pvalCutoff$Log2Rat))))),]
#         belowLODR <- subset(pval,
#                             (Log2Rat >= logFCcut) & 
#                               (MnCnt <= estLODRraw ) & 
#                               (pvals <= p.thresh)) 
      }
      
      cat(paste("\nFlagged",as.character(dim(belowLODR)[1]),
                "DE Features that may be false negatives\n",
                "Pvalues are significant for chosen FDR, but fold change\n",
                "measurements are less than the fold change threshold\n"))
      ## Return portion of transcript list that is below the LODR for the chosen FC and p.thresh
      nam <- paste("belowLODR",as.character(FC[i]),sep = ".")
      expDat$Results$belowLODR <- belowLODR
      names(expDat$Results)[which(names(expDat$Results) == "belowLODR")] <- nam
    }else{
      cat(paste("\n LODR estimate is undefined for abs(log2(fold change)) =",
                round(logFCcut,digits=3),"(Ratio = ",Ratio[i],")\n"))
    }
    
  }
  #expDat$ResultspvalCutoff <- pvalCutoff
  return(expDat)
    
}