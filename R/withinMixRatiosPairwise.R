withinMixRatios <- function(exDat){
  # calculate deviation of relative abundances of transcripts within a mix from
  #the nominal concentration
  
 
  colScale <- exDat$plotInfo$colScale
  fillScale <- exDat$plotInfo$fillScale
  datType <- exDat$sampleInfo$datType
  sampleNames <- exDat$sampleNames
  
  playDat <- merge(exDat$normERCCDat, exDat$erccInfo$idColsSRM)
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
  for (j in 1:nlevels(playDat$Ratio)){
    subDyn <- subset(playDat, 
                     (Ratio == levels(playDat$Ratio)[j]))
    #print(nrow(subDyn))
    nominalMix1 <- NULL
    nominalMix2 <- NULL
    observedMix1 <- NULL
    observedMix2 <- NULL
    aveCtl1<-NULL
    aveCtl2<-NULL
    FeaturePair <- NULL
    Ratio <- NULL
    GCrat <- NULL
    Lengthrat <- NULL
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
      GCrat[k-1] <- subDyn$GC[k-1]/subDyn$GC[k]
      Lengthrat[k-1] <- subDyn$Length[k-1]/subDyn$Length[k]
    }
    
    
    addRes1 <- data.frame(Feature = FeaturePair, Ratio = Ratio,
                          Nominal = nominalMix1, Observed = observedMix1, 
                          AveCtl = aveCtl1,
                          Sample = sampleNames[1], GCrat = GCrat, 
                          Lengthrat = Lengthrat)
    addRes2 <- data.frame(Feature = FeaturePair, Ratio = Ratio,
                          Nominal = nominalMix2, Observed = observedMix2,
                          AveCtl = aveCtl2,
                          Sample = sampleNames[2], GCrat = GCrat, 
                          Lengthrat = Lengthrat)
    #addRes <- cbind(subDyn, addRes)
    addRes <- rbind(addRes1, addRes2)
    perERCCRatios <- rbind(perERCCRatios, addRes)
    #print(perERCCRatios)
  }
  
  perERCCRatios$NomFactor <- as.factor(round(abs(perERCCRatios$Nominal)))
  
  betwCtlRats <- ggplot(perERCCRatios, aes(x = NomFactor, y = Observed)) + 
    geom_boxplot(aes(x= NomFactor), alpha = 0.5) + 
    geom_point(aes(colour = Ratio,shape = Sample),size = 5,alpha = 0.7) + 
    ylab("Observed Log2 Ratios for pairs of ERCCs") + 
    xlab("Nominal Log2 Ratios for pairs of ERCCs") + 
    facet_grid(~Sample) + colScale + fillScale + theme_bw() #+ geom_text(aes(label = gsub("ERCC-00", "",Feature)),angle = 315)
  
  residRats <- ggplot(perERCCRatios, aes(x = NomFactor, y = Observed-Nominal)) +
    geom_boxplot(aes(x= NomFactor), alpha = 0.5) + 
    geom_point(aes(colour = Ratio,shape = Sample),size = 5,alpha = 0.7) + 
    ylab(paste0("Residuals (Observed - Nominal) of Log2 Ratios \n", 
         "for Pairs of ERCC Controls Within Each Mixture")) + 
    xlab(paste0("Nominal Log2 Ratios \n",
                "for Pairs of ERCC Controls Within Each Mixture")) +
    facet_grid(~Sample) + 
    colScale + fillScale + geom_abline(aes(intercept = 0, slope =0)) + theme_bw()
  
  perERCCRatios$Resid <- (perERCCRatios$Observed - perERCCRatios$Nominal)
  getVertices = by(perERCCRatios, perERCCRatios[,"Ratio"], function(myHull) chull(myHull$AveCtl,myHull$Resid))
  for (t in 1:(length(getVertices))){
    findVerts = perERCCRatios[which(perERCCRatios$Ratio == names(getVertices[t])),]
    # handle the 1st time through
    if (t==1){
      getVerts = findVerts[getVertices[[t]],]
    }else{	
      getVerts = rbind(getVerts, findVerts[getVertices[[t]],] )		
    }	
  }

  medList <- as.data.frame(tapply(perERCCRatios$Resid,perERCCRatios[,"Ratio"],median))
  medList$Ratio <- row.names(medList)
  names(medList)[1]<- "medianResid"
  medList$Ratio <- as.factor(medList$Ratio)
  medList$medianResid <- as.numeric(medList$medianResid)
  
  ratVerts <- ggplot()+
    geom_polygon(data = getVerts, aes(x = AveCtl, 
                                      y = Resid, 
                                      fill = Ratio), alpha = 0.15) +
    geom_point(data = perERCCRatios, aes(x = AveCtl,
                                         y = Resid,
                                         colour = Ratio, 
                                         shape = Sample), 
               size = 5, alpha = 0.5) + 
    ylab(paste0("Residuals (Observed - Nominal) of Log2 Ratios \n",
                "for Pairs of ERCC Controls Within Each Mixture")) + 
    xlab(paste0("Average Log2 Abundance \n",
                "for Pairs of ERCC Controls Within Each Mixture")) +
    geom_hline(aes(yintercept = 0), colour = "grey50") +
    geom_hline(data = medList, aes(yintercept = medianResid, colour = Ratio)) +
    colScale + fillScale + theme_bw()
  
  
  violinrats <- ggplot(perERCCRatios, aes(x = NomFactor, y = Observed-Nominal))+
    geom_violin(aes()) + 
    ylab("Residuals for ERCC pairwise ratios (Log2(Observed) - Log2(Nominal))")+
    xlab("Nominal Log2 Ratios for pairs of ERCCs") + 
    facet_grid(~Sample) + geom_abline(aes(slope = 0, intercept = 0)) + theme_bw()
    
  
  #multiplot(betwCtlRats,residCtlRats,cols=1)
  
  #ggplot(perERCCRatios)+ geom_point(aes(x = Observed, y = Nominal))
  exDat$Results$perERCCRatios <- perERCCRatios
  exDat$Figures$betwCtlRats <- betwCtlRats
  exDat$Figures$residRats <- residRats
  exDat$Figures$violinrats <- violinrats
  exDat$Figures$ratVerts <- ratVerts
  return(exDat)
  
}