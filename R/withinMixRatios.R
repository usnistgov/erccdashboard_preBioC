withinMixRatios <- function(expDat){
  # calculate deviation of relative abundances of transcripts within a mix from
  #the nominal concentration
  
  dynRangeDat <- expDat$Results$dynRangeDat
  colScale <- expDat$plotInfo$colScale
  fillScale <- expDat$plotInfo$fillScale
  
  dynRangeDat_o <- dynRangeDat[order(dynRangeDat$Ratio,
                                     dynRangeDat$Sample,
                                     dynRangeDat$Conc,decreasing = T),]
  perERCCRatios <- NULL
  nominal <- NULL
  observed <- NULL
  Ratio <- NULL
  Sample <- NULL
  Feature <- NULL
  for (i in 1:nlevels(dynRangeDat_o$Sample)){
    #subDynSample <- subset(dynRangeDat_o, Sample == 
    #                   as.character(levels(dynRangeDat_o$Sample))[i])
    for (j in 1:nlevels(dynRangeDat_o$Ratio)){
      subDyn <- subset(dynRangeDat_o, 
                       (Ratio == levels(dynRangeDat_o$Ratio)[j]) &
                         (Sample == levels(dynRangeDat_o$Sample)[i]))
      print(nrow(subDyn))
      nominal <- NULL
      observed <- NULL
      for (k in 1:nrow(subDyn)){
        #print(nrow(subDyn))
        
        nominal[k] <- subDyn$Conc[1] - subDyn$Conc[k]
        observed[k] <- subDyn$value.Ave[1] - subDyn$value.Ave[k]
        
      }
      Feature <- as.character(subDyn$Feature)
      Ratio <- as.character(subDyn$Ratio)
      Sample <- as.character(subDyn$Sample)
      print(nominal)
      print(observed)
      addRes <- data.frame(Feature = Feature, Ratio = Ratio, Sample = Sample,
                           Nominal = nominal, Observed = observed)
      #addRes <- cbind(subDyn, addRes)
      perERCCRatios <- rbind(perERCCRatios, addRes)
      print(perERCCRatios)
    }
  }
  
  ggplot(perERCCRatios) + geom_point(aes(x = Nominal, y = Observed-Nominal, colour = Ratio, shape = Sample),size = 5)

  ggplot(perERCCRatios, aes(x = Nominal, y = Observed-Nominal)) + geom_boxplot(aes()) + geom_point(aes(colour = Ratio, shape = Sample),size = 5) + ylab("Observed Ratio - Nominal Ratio for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(~Ratio)
  
  ggplot(perERCCRatios, aes(x = Nominal, y = Observed)) + geom_point(aes(colour = Ratio, shape = Sample),size = 5) + ylab("Observed Ratios for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(~Ratio)
  
  
  betwCtlRats <- ggplot(perERCCRatios, aes(x = Nominal, y = Observed)) + geom_point(aes(colour = Ratio, shape = Sample),size = 5,alpha = 0.7) + ylab("Observed Ratios for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(Sample~Ratio) + geom_text(aes(label = gsub("ERCC-00", "",Feature)),hjust = -0.5) + geom_abline(aes(slope = 1, intercept = 0)) + colScale + fillScale
  
  residCtlRats <- ggplot(perERCCRatios, aes(x = Nominal, y = Observed-Nominal)) + geom_boxplot(aes(fill = Ratio), alpha = 0.5) + geom_point(aes(colour = Ratio, shape = Sample),size = 5,alpha = 0.7) + ylab("Observed Ratios for pairs of ERCCs") + xlab("Nominal Ratios for pairs of ERCCs") + facet_grid(Sample~Ratio) + geom_text(aes(label = gsub("ERCC-00", "",Feature)),angle = 315) + geom_abline(aes(slope = 0, intercept = 0))+ colScale + fillScale
  
  multiplot(betwCtlRats,residCtlRats,cols=1)
  
  expDat$Figures$betwCtlRats <- betwCtlRats
  expDat$Figures$residCtlRats <- residCtlRats
    
}