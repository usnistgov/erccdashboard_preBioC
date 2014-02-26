arrayAnalysis<- function(expDat){
  
  sampleInfo <- expDat$sampleInfo
  erccInfo <- expDat$erccInfo
  plotInfo <- expDat$plotInfo
  
  library("limma")
  COHarray <- read.table("../erccdashboard201309/data/QCmicroarrays_COH_UTSW/COHdata110909.txt",header = T)
  COHarrayAB <- COHarray[-c(11:16)]
  y <- COHarrayAB[c(1,2,5:10)]
  y$PROBE_ID<- as.character(y$PROBE_ID)
  idxERCC <- grep("ERCC-",y$TargetID)
  
  # copy the ERCC target ids over to the probe ID column
  y$PROBE_ID[idxERCC] <- as.character(y$TargetID[idxERCC])
  
  row.names(y)<- y$PROBE_ID
  y <- y[-c(1,2)]
  
  ylog <- log2(y)
  design <- cbind(Grp1=1,Grp1vs2=c(1,1,1,0,0,0))
  
  fit <- lmFit(ylog,design)
  fit <- eBayes(fit)
  res <- topTable(fit,number = dim(ylog)[1],coef = 2)
  ERCCres<- res[grep("ERCC-",res$ID),]
  
  plot(ERCCres$AveExpr,ERCCres$logFC)
  
  pval.res <- data.frame( Feature = ERCCres$ID, MnCnt = ERCCres$AveExpr, Pval = ERCCres$adj.P.Val, Fold = )
 
  write.csv(pval.res,file=paste(sampleInfo$filenameRoot,"ERCC Pvals.csv"),
            row.names = F)
}