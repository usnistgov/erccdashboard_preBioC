testDEArray <- function(expDat){
  library(limma)
  y <- expDat$Transcripts
  choseFDR <- expDat$sampleInfo$choseFDR
  if(is.null(choseFDR)){
    p.thresh <- expDat$Results$p.thresh
  }
  
  #normy <- sweep(y[-c(1)], 2, expDat$libeSize/1000,"/")
  #y <- cbind(y[c(1)],normy)
  erccInfo <- expDat$erccInfo
  sampleInfo <- expDat$sampleInfo
  row.names(y) <- y$Feature
  y <- y[-c(1)]
  
  
  ylog <- log2(y)
  
  wts <- log(expDat$libeSize)
  
  design <- cbind(Grp1=1,Grp1vs2=c(1,1,1,0,0,0))
  
  fit <- lmFit(ylog,design,weights=wts)
  ###fit <- lmFit(ylog,design)
  
  fit <- eBayes(fit)
  
  res <- topTable(fit,sort.by="none",number = dim(ylog)[1],coef = 2)
  #print(head(res))
  ### generate qvals
  if(!is.null(choseFDR)){
    
    pval <- res$adj.P.Val
    qobj <- qvalue(pval)
    res$qvals <- qobj$qvalues
    
    if(any(res$qvals<choseFDR)){
      p.thresh<-max(res$adj.P.Val[res$qvals<choseFDR])
    }
  }
  res$Feature <- row.names(res)
  
  ERCCres<- res[grep("ERCC-",row.names(res)),]
  
  #plot(ERCCres$AveExpr,ERCCres$logFC)
  
  erccFC <- data.frame(Feature = erccInfo$idColsSRM$Feature, 
                       FC = round(erccInfo$idColsSRM$Conc1/
                                    erccInfo$idColsSRM$Conc2,digits=3))
  
  pval.res <- data.frame( Feature = row.names(ERCCres), MnCnt = 2^(ERCCres$AveExpr), 
                          Pval = ERCCres$adj.P.Val,
                          Fold = erccFC$FC[match(row.names(ERCCres), erccFC$Feature,
                                                 nomatch=0)])
  write.csv(pval.res, paste(sampleInfo$filenameRoot, "ERCC Pvals.csv"),
            row.names = F)
  
  expDat$Results$limma.res <- res
  expDat$Results$ERCC.pval <- pval.res
  expDat$Results$p.thresh <- p.thresh
  return(expDat)
}