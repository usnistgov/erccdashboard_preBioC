meltExpDat <- function(expDat, cnt, designMat){
  sampleInfo <- expDat$sampleInfo
  libeSize <- expDat$sampleLibeSums/(10^6)
  datNames <- colnames(designMat)[-1]
  sample1 <- expDat$sampleNames[1]
  sample2 <- expDat$sampleNames[2]
  
  datCols = cnt[-c(1)]
  libAdjust = sweep(datCols, 2, libeSize,"/")
  sampleLibeDataNorm = cbind(cnt[c(1)],libAdjust)
  myDataERCC = sampleLibeDataNorm
  expressDat = myDataERCC[-c(1)] 
  sampleNameList = c(sample1,sample2)
  # added 131224
  dat = cbind(Feature = myDataERCC[c(1)], Ratio = "Endo",
                               myDataERCC[-c(1)])
  
  dat$Feature <- as.factor(as.character(
    dat$Feature))
  dat$Ratio <- as.character(dat$Ratio)
  idx1 <- suppressWarnings(which(as.character(dat$Feature) == 
                 as.character(expDat$idColsAdj$Feature)))
  dat$Ratio[idx1] <- as.character(expDat$idColsAdj[idx1,c(4)])
  dat$Ratio <- as.factor(dat$Ratio)
  #dat = merge(expDat$idColsAdj[c(1,4)],myDataERCC)
  
  #dat$Feature <- as.factor(as.character(
  #  dat$Feature))
  #dat$Ratio <- factor(as.character(dat$Ratio),
  #                                     levels = as.character(
  #                                       expDat$sampleInfo$FCcode$Ratio))
  head(dat)
  #expressDat_l <- melt(dat)
  
  colAdd <- colsplit(expressDat_l$variable,pattern="_",names=datNames)
  colAdd = as.data.frame(lapply(colAdd,as.factor))
  expressDat_l<- expressDat_l[,-c(3)]
  colnames(expressDat_l)[3]<-"NormCounts"
  expressDat_l <- cbind(expressDat_l, colAdd)
  
  #expressDat_l$Sample <- as.factor(as.character(expressDat_l$Sample))
  #expressDat_l$Replicate <- as.factor(as.character(expressDat_l$Replicate))
  expDat$expressDat_l <- expressDat_l
  return(expDat)
}