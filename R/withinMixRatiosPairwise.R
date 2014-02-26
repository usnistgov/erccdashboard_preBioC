withinMixRatios <- function(expDat){
  # calculate deviation of relative abundances of transcripts within a mix from
  #the nominal concentration
  
  #dynRangeDat <- expDat$Results$dynRangeDat
  colScale <- expDat$plotInfo$colScale
  fillScale <- expDat$plotInfo$fillScale
  sampleNames <- expDat$sampleNames
  #dynRangeDat_o <- dynRangeDat[order(dynRangeDat$Ratio,
  #                                   dynRangeDat$Sample,
  #                                   dynRangeDat$Conc,decreasing = T),]
  playDat <- merge(expDat$normERCCDat, expDat$erccInfo$idColsSRM)
  playDat$Feature <- as.factor(as.character(playDat$Feature))
  
  playDat$Mean1 <- apply(playDat[grep(sampleNames[1],colnames(playDat))],1,mean)
  playDat$Mean2 <- apply(playDat[grep(sampleNames[2],colnames(playDat))],1,mean)
  
  playDat <- playDat[order(playDat$Ratio,playDat$Conc1),]
  playDat <- playDat[-grep(paste0(sampleNames[1],"|",sampleNames[2]),colnames(playDat))]
  #playDat_l <- melt(playDat[-grep(paste0(sampleNames[1],"|",sampleNames[2]),colnames(playDat))])
  perERCCRatios <- NULL
  nominal <- NULL
  observed <- NULL
  Ratio <- NULL
  Sample <- NULL
  Feature <- NULL
  #for (i in 1:2)){
    #subDynSample <- subset(dynRangeDat_o, Sample == 
    #                   as.character(levels(dynRangeDat_o$Sample))[i])
  for (j in 1:nlevels(playDat$Ratio)){
    subDyn <- subset(playDat, 
                     (Ratio == levels(playDat$Ratio)[j]))
    print(nrow(subDyn))
    nominalMix1 <- NULL
    nominalMix2 <- NULL
    observedMix1 <- NULL
    observedMix2 <- NULL
    aveCtl1<-NULL
    aveCtl2<-NULL
    FeaturePair <- NULL
    Ratio <- NULL
    #Sample <- NULL
    for (k in 2:nrow(subDyn)){
      #print(nrow(subDyn))
      nominalMix1[k-1] <- log2(subDyn$Conc1[k-1]) - log2(subDyn$Conc1[k])
      nominalMix2[k-1] <- log2(subDyn$Conc2[k-1]) - log2(subDyn$Conc2[k])
      aveCtl1[k-1]<- log2((subDyn$Mean1[k-1] + subDyn$Mean1[k])/2) 
      aveCtl2[k-1]<- log2((subDyn$Mean1[k-1] + subDyn$Mean1[k])/2) 
      observedMix1[k-1] <- log2(subDyn$Mean1[k-1]) - log2(subDyn$Mean1[k])
      observedMix2[k-1] <- log2(subDyn$Mean2[k-1]) - log2(subDyn$Mean2[k])
      FeaturePair[k-1] <- paste(subDyn$Feature[k],
                                subDyn$Feature[k-1],sep = ".")
      Ratio[k-1] <- as.character(subDyn$Ratio[k-1])
      #Sample[k-1]<- as.character(subDyn$Sample[k-1])
      
    }
    #Feature <- as.character(subDyn$Feature)
    #Ratio <- as.character(subDyn$Ratio)
    #Sample <- as.character(subDyn$Sample)
    
    addRes1 <- data.frame(Feature = FeaturePair, Ratio = Ratio,
                          Nominal = nominalMix1, Observed = observedMix1, 
                          AveCtl = aveCtl1,
                          Sample = sampleNames[1])
    addRes2 <- data.frame(Feature = FeaturePair, Ratio = Ratio,
                          Nominal = nominalMix2, Observed = observedMix2,
                          AveCtl = aveCtl2,
                          Sample = sampleNames[2])
    #addRes <- cbind(subDyn, addRes)
    addRes <- rbind(addRes1, addRes2)
    perERCCRatios <- rbind(perERCCRatios, addRes)
    print(perERCCRatios)
  }
  #}
  perERCCRatios$NomFactor <- as.factor(round(perERCCRatios$Nominal))
  
  
  betwCtlRats <- ggplot(perERCCRatios, aes(x = NomFactor, y = Observed)) + 
    geom_boxplot(aes(x= NomFactor), alpha = 0.5) + 
    geom_point(aes(colour = Ratio,shape = Sample),size = 5,alpha = 0.7) + 
    ylab("Observed Log2 Ratios for pairs of ERCCs") + 
    xlab("Nominal Log2 Ratios for pairs of ERCCs") + 
    facet_grid(~Sample) + colScale + fillScale #+ geom_text(aes(label = gsub("ERCC-00", "",Feature)),angle = 315)
  
  print(betwCtlRats)
  
  residRats <- ggplot(perERCCRatios, aes(x = NomFactor, y = Observed-Nominal)) +
    geom_boxplot(aes(x= NomFactor), alpha = 0.5) + 
    geom_point(aes(colour = Ratio,shape = Sample),size = 5,alpha = 0.7) + 
    ylab("Observed/Nominal Log2 Ratios for pairs of ERCCs") + 
    xlab("Nominal Log2 Ratios for pairs of ERCCs") + facet_grid(~Sample) + 
    colScale + fillScale + geom_abline(aes(intercept = 0, slope =0))#+ geom_text(aes(label = gsub("ERCC-00", "",Feature)),angle = 315)
  
  print(residRats)
  
  ERCCEffects <- ggplot(dynRangeDat, aes(x=Conc, y = (value.Ave-Conc)-median(value.Ave-Conc)))+ geom_boxplot()+geom_point(aes(colour = Ratio, shape = Sample),size = 5, alpha = 0.5) + colScale + fillScale
  
  ggplot(perERCCRatios) + geom_point(aes(x = Nominal, y = Observed-Nominal, colour = Ratio, shape = Sample),size = 5)

  ggplot(perERCCRatios, aes(x = Nominal, y = Observed-Nominal)) + geom_boxplot(aes()) + geom_point(aes(colour = Ratio, shape = Sample),size = 5) + ylab("Observed Ratio - Nominal Ratio for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(~Ratio)
  
  ggplot(perERCCRatios, aes(x = Nominal, y = Observed)) + geom_point(aes(colour = Ratio, shape = Sample),size = 5) + ylab("Observed Ratios for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(~Ratio)
  
  
  betwCtlRats <- ggplot(perERCCRatios, aes(x = Nominal, y = Observed)) + geom_point(aes(colour = Ratio, shape = Sample),size = 5,alpha = 0.7) + ylab("Observed Ratios for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(Sample~Ratio) + geom_text(aes(label = gsub("ERCC-00", "",Feature)),hjust = -0.5) + geom_abline(aes(slope = 1, intercept = 0)) + colScale + fillScale
  
  residCtlRats <- ggplot(perERCCRatios, aes(x = Nominal, y = Observed-Nominal)) + geom_boxplot(aes(fill = Ratio), alpha = 0.5) + geom_point(aes(colour = Ratio, shape = Sample),size = 5,alpha = 0.7) + ylab("Observed Ratios for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(Sample~Ratio) + geom_text(aes(label = gsub("ERCC-00", "",Feature)),angle = 315) + geom_abline(aes(slope = 0, intercept = 0))+ colScale + fillScale
  
  
  multiplot(betwCtlRats,residCtlRats,cols=1)
  
  ggplot(perERCCRatios)+ geom_point(aes(x = Observed, y = Nominal))
  
  expDat$Figures$betwCtlRats <- betwCtlRats
  expDat$Figures$residCtlRats <- residCtlRats
    
}