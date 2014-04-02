# Plots with target ratios and R adjusted ratios, MA plots and Ratio Summaries
maConcPlot <-function(expDat, LODR.annot.ERCC, alphaPoint = 0.8,
                      r_mAdjust = T, replicate = T){
  
  ReplicateName = "Rep"
  # Melt the data
  expDat <- meltExpDat(expDat, cnt = expDat$Transcripts, 
                       designMat = expDat$designMat)
  countPair <- expDat$expressDat_l
  
  sampleInfo <- expDat$sampleInfo
  myYLim = sampleInfo$myYLimMA
  myXLim = sampleInfo$myXLim
  siteName = sampleInfo$siteName
  analysis = sampleInfo$analysis
  filenameRoot = sampleInfo$filenameRoot
  sample1 = expDat$sampleNames[1]
  sample2 = expDat$sampleNames[2]
 
  
  idCols = expDat$idColsAdj
  r_m.res = expDat$r_m.res
  
  cutoffs = LODR.annot.ERCC$cutoffs
  FCcode = sampleInfo$FCcode
  legendLabels=sampleInfo$legendLabels
  avexlabel = expDat$ERCCxlabelAve
  spikeFraction = expDat$spikeFraction
  
  r_m.mn = exp(r_m.res$r_m.mn)
  r_m.UL = exp(r_m.res$r_m.upper)
  r_m.LL = exp(r_m.res$r_m.lower)
  
  theme_update(legend.justification=c(1,0), legend.position=c(1,0))
  
  plotParams <- plotAdjust(expDat)
  colScale <- plotParams$colScale
  fillScale <- plotParams$fillScale
  
  library(gridExtra)
  
  if (length(cutoffs) > 0){
    cat("\nLODR estimates are available to code ratio-abundance plot\n")  
  }
  
  counts1 = countPair$NormCounts[which(countPair$Sample == sample1)]
  counts2 = countPair$NormCounts[which(countPair$Sample == sample2)]
  
  if (length(counts1) != length(counts2)){
    stop("Uneven replication for the sample pair")
  }
 
  
  if (replicate == T){
    Mdat = log2(counts1) - log2(counts2)
    Adat = log2((counts1 + counts2)/2)
    colIdx <- which(colnames(countPair) == ReplicateName)
    maData = data.frame(Feature = countPair$Feature[which(countPair$Sample == 
                                                            sample1)], 
                        Ratio = countPair$Ratio[which(countPair$Sample ==
                                                        sample1)],
                        M = Mdat, A = Adat, Replicate = 
                          countPair[which(countPair$Sample == sample1),colIdx])
  }
  else{
    stop("Replicate samples are needed")
  }
    
  names(maData)[3:4] = c("M","A")
  expDat$maData <- maData

  ### Now subset maData and continue with just ERCCs
  maData <- subset(maData, subset = Ratio != "Endo")
  maData$Ratio <- factor(as.character(maData$Ratio),levels=FCcode$Ratio)
  maData$Feature <- as.factor(as.character(maData$Feature))
  #countPair <- subset(expDat$expressDat_l, subset = Ratio != "Endo")
  #countPair$Ratio <- as.factor(as.character(countPair$Ratio))
  
  maData$Nominal = FCcode$FC[1]
  for (i in 2:nlevels(FCcode$Ratio)){
    maData$Nominal[which(maData$Ratio == FCcode$Ratio[i])] = FCcode$FC[i]
  }
 
    
  if(r_mAdjust == T){
    maData$Empirical = maData$Nominal/r_m.mn  
  }
  xlabel = xlab(avexlabel)
  #if(replicate == T){
   # Use Average Concentration values (instead of Average Signals)
    FeatureA = data.frame(Feature= idCols$Feature, A = log2((idCols$Conc1 +
                                                               idCols$Conc2)/2))
   for (i in 1:nlevels(FeatureA$Feature)){
     maData$A[which(maData$Feature == as.character(FeatureA$Feature[i]))] =
       FeatureA$A[i]  
   }
    maDataAve <- data.frame(tapply(maData$M,maData[,"Feature"],mean))
    
    maDataAveSD = cbind(row.names(maDataAve),
                      maData$Ratio[match(row.names(maDataAve),
                                         table=maData$Feature)],
                      maData$Empirical[match(row.names(maDataAve),
                                             table=maData$Feature)],
                      maData$A[match(row.names(maDataAve),
                                     table=maData$Feature)], maDataAve[c(1)])
  
  names(maDataAveSD)[1:5] = c("Feature","Ratio","Empirical","A","M.Ave")
  
  maDataSD <- as.vector(tapply(maData$M,maData[,"Feature"],sd))
  maDataAveSD$M.SD = maDataSD
  
  cutERCCs = unique(maDataAveSD$Feature[which(is.na(maDataAveSD$M.SD))])
    
  if(length(cutERCCs) != 0){
    cat(paste("\nThese ERCCs were not included in the ratio-abundance plot, \n",
              "because not enough non-zero replicate measurements of these \n",
              "controls were obtained for both samples:\n"))
    
    cat(paste(as.character(cutERCCs),collapse="\n"))
    
    maDataAveSD = maDataAveSD[-(which(maDataAveSD$Feature %in% cutERCCs)),]
    searchCut = paste(cutERCCs,collapse="|")
    maData = maData[-(grep(searchCut,maData$Feature)),]
  }
    
    maDataAveSD$Feature = as.factor(as.character(maDataAveSD$Feature))  
    maData$Feature = as.factor(as.character(maData$Feature))
    #write.csv(maData,file = paste(filenameRoot,"maDataFinite.csv",sep = "."))
    
    # Estimate SD for all ratio measurements
    # For the ERCCs that are plotted take the log ratio data and subtract the log 
    # nominal ratio, this will normalize the data
    normRats = maData$M - log2(maData$Nominal)
    #print(normRats)
    
    # Take the sd of all of the normalized log ratio data to find a global SD for the ERCCs at this site
    sdGlobal = sd(normRats)
    #cat("\nGlobal Ratio SD for this sample pair is: ")
    #cat(sdGlobal, "\n")
    
    ratVarDat <- maDataAveSD
    
    ratVarDat$A <- log2((2^(ratVarDat$A))/(spikeFraction))
    
    
    
    
    sdRatioplot <- ggplot() + geom_point(data = ratVarDat, aes(x = A, y = M.SD, 
                                                               colour = Ratio),
                                         size = 5, alpha = alphaPoint)+colScale 
      
    
    lmfit <- lm(formula = M.SD ~ 1, data = ratVarDat, weights = 1/((M.SD)^2))
    
    stdevCoef <- summary(lmfit)$coefficients
    
    params<-list(Mo = 0.1, lambda = 2) # setup default params.
    f = coef(lmfit)[[1]]  
    #print("Using default parameters")
    #print(params)
    #print(f)
    
    while(TRUE){
      
      expmod<-NULL
      # does not stop in the case of error
      try(expmod<-nls(formula=M.SD ~ Mo*exp(-lambda * A) + f, data = ratVarDat,
                      start = params), silent = TRUE)
      if(!is.null(expmod))break; # if nls works, then quit from the loop
      #print("Try changing initial guess for parameters")
      params<-list(Mo = sdGlobal, lambda = 0.2)# change the params for nls
      #print(params)
      
    }
    
    if(!is.null(expmod)){
      stdevCoef <- rbind(stdevCoef, summary(expmod)$coefficients)
    } 
    
    #print(stdevCoef)
    row.names(stdevCoef)<- c("Minimum SD Estimate (Vmin)", 
                             "Maximum SD Estimate (Vmax)",
                             "Lambda")
    print(stdevCoef[,c(1:2)])
    
    
    # add fitted curve
    lineDat <- data.frame("xdat" = ratVarDat$A,"ydat" = 
                            predict(lmfit, x = ratVarDat$A)) 
    if(!is.null(expmod)){
      fitDat <- data.frame("xdat" = ratVarDat$A,"ydat" = predict(expmod,
                                                                 list(x = ratVarDat$A))) 
    } 
    sdRatioplotFit <- sdRatioplot + geom_line(data = lineDat,aes(x=xdat,y=ydat))+
      ylab("sd(Log2 Ratios of Counts)") + xlabel + 
      theme(legend.justification=c(1,1), legend.position=c(1,1))
    
    
    if(!is.null(expmod)){
      sdRatioplotFit <- sdRatioplotFit + 
        geom_line(data = fitDat, aes(x = xdat, y = ydat)) + xlabel + 
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      
    } 
    
    
    if (length(cutoffs)>0){
      maDataAveSD$LODR = "below"
      FCcodeC = FCcode[-c(which(FCcode$FC == 1)),]
      for (i in 1:length(cutoffs)){
        maDataAveSD$LODR[which((maDataAveSD$A > cutoffs[i])&
                                 (maDataAveSD$Ratio ==
                                    FCcodeC$Ratio[i]))] = "above"
      }
      
      maDataAveSD$LODR <- as.factor(maDataAveSD$LODR)
      rm_dat = data.frame("Mean r_m" = signif(r_m.mn,digits = 4),
                          "Lower Limit" = signif(r_m.LL,digits = 4), 
                          "Upper Limit" = signif(r_m.UL,digits = 4))
      colnames(rm_dat) <- c("Mean", expression(paste("95% CI Lower Bound")), 
                            expression(paste("95% CI Upper Bound")))
      rownames(rm_dat) <- expression(r[m])
      
      # Plot ratio-signal data coding points below the LODR with open circles
      maPlot <- ggplot(maDataAveSD, aes(x = A, y = M.Ave) ) + 
        geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD, 
                          colour = Ratio),alpha = 0.3) + 
        geom_point(aes(colour = Ratio),size = 5, alpha = alphaPoint) +
        geom_point(data = subset(maDataAveSD, (LODR == "below")),
                   colour = "white",size = 2.5) + 
        geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                   linetype = "longdash",alpha = 0.7) + 
        ylab("Log2 Ratio of Counts") + xlabel + 
        coord_cartesian(xlim = myXLim, ylim = myYLim) + colScale + 
        annotation_custom(tableGrob(rm_dat,parse=T, 
                                    gpar.corefill = gpar(fill = "grey85",
                                                        col = "white"), 
                                    gpar.rowfill = gpar(fill = "grey80",
                                                        col = "white"),
                                    gpar.colfill = gpar(fill = "grey80",
                                                        col = "white")), xmin = 
                            quantile(maDataAveSD$A,probs=0.25), 
                          xmax = max(maDataAveSD$A),
                          ymin = (myYLim[2]) - 0.25*myYLim[2], 
                          ymax = myYLim[2]) + theme(legend.justification=c(1,0),
                                                    legend.position=c(1,0)) 
    }else{
      # Plot ratio-signal data without LODR coding
      maPlot <- ggplot(maDataAveSD, aes(x = A, y = M.Ave, colour = Ratio) ) + 
        geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD), 
                      alpha = 0.3) + geom_point(size = 5, alpha = alphaPoint) + 
        geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                   linetype = "longdash",alpha = 0.7) + 
        ylab("Log2 Ratio of Counts") + xlabel + colScale+ 
        coord_cartesian(xlim = myXLim, ylim = myYLim) + 
        theme(legend.justification=c(1,0), legend.position=c(1,0)) 
    }
  
  expDat$Figures$plotRatioAnnot <- maPlot
  expDat$sdGlobal <- sdGlobal
  expDat$modRatVar <- stdevCoef
  expDat$Figures$plotSDratio <- sdRatioplotFit
  expDat$maDataAveSD
  return(expDat)
}