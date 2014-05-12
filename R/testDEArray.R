testDEArray <- function(exDat){
  #library(limma)
  y <- exDat$Transcripts
  choseFDR <- exDat$sampleInfo$choseFDR
  if(is.null(choseFDR)){
    p.thresh <- exDat$Results$p.thresh
  }
  
  #normy <- sweep(y[-c(1)], 2, exDat$normFactor/1000,"/")
  #y <- cbind(y[c(1)],normy)
  erccInfo <- exDat$erccInfo
  sampleInfo <- exDat$sampleInfo
  row.names(y) <- y$Feature
  y <- y[-c(1)]
  
  
  ylog <- log2(y)
  if(is.null(exDat$normFactor)&(exDat$sampleInfo$isNorm==TRUE)){
    wts = NULL
    cat("\nisNorm is TRUE, array data is already normalized\n")
  }else{
    wts <- log(exDat$normFactor)  
  }
  
  if(odd(ncol(y))) stop("\nUneven number of replicates for the two sample types\n")
  
  design <- cbind(Grp1=1,Grp1vs2=c(rep(x=1,times=ncol(y)/2), rep(x=0,times=ncol(y)/2)))
  
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
  
  pval.res <- data.frame( Feature = row.names(ERCCres), MnSignal = 2^(ERCCres$AveExpr), 
                          Pval = ERCCres$adj.P.Val,
                          Fold = erccFC$FC[match(row.names(ERCCres), erccFC$Feature,
                                                 nomatch=0)])
  write.csv(pval.res, paste(sampleInfo$filenameRoot, "ERCC Pvals.csv"),
            row.names = F)
  
  exDat$Results$limma.res <- res
  exDat$Results$ERCC.pval <- pval.res
  exDat$Results$p.thresh <- p.thresh
  return(exDat)
}