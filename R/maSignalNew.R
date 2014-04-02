#' Generate MA plots with or without annotation using LODR estimates 
#'
#' @param expDat      list, contains input data and stores analysis results
#' @param alphaPoint  numeric value, for alpha (transparency) for plotted points,
#'                    range is 0 - 1
#' @param r_mAdjust   default is TRUE, if FALSE then the r_m estimate will not
#'                    used to offset dashed lines for empirical ratios on figure
#' @param replicate   default is TRUE, if FALSE then error bars will not be
#'                    produced
#'                    
#' 
#' @export

# Plots with target ratios and R adjusted ratios, MA plots and Ratio Summaries
maSignal <-function(expDat, alphaPoint = 0.8, r_mAdjust = T, replicate = T){
  
  #ReplicateName = "Rep"
  # Melt the data
  #expDat <- meltExpDat(expDat, cnt = expDat$Transcripts, 
  #                    designMat = expDat$designMat)
  #countPair <- expDat$expressDat_l
  sampleInfo <- expDat$sampleInfo
  erccInfo <- expDat$erccInfo
  plotInfo <- expDat$plotInfo
  
  cnt <- expDat$Transcripts
  designMat <- expDat$designMat
  sampleInfo <- expDat$sampleInfo
  libeSize <- expDat$libeSize
  datNames <- colnames(designMat)[-1]
  sample1 <- expDat$sampleNames[1]
  sample2 <- expDat$sampleNames[2]
  
  datCols = cnt[-c(1)]
  libAdjust = sweep(datCols, 2, libeSize,"/")
  sampleLibeDataNorm = cbind(cnt[c(1)],libAdjust)
  myDataERCC = sampleLibeDataNorm
  expressDat = myDataERCC[-c(1)] 
  sampleNameList = c(sample1,sample2)
  
  dat = cbind(Feature = myDataERCC[c(1)], Ratio = "Endo",
              myDataERCC[-c(1)])
  
  dat$Feature <- as.factor(as.character(
    dat$Feature))
  dat$Ratio <- as.character(dat$Ratio)
  #idx1 <- suppressWarnings(which(as.character(dat$Feature) == 
  #                                 as.character(expDat$idColsAdj$Feature)))
  idxERCC <- grep("ERCC-",dat$Feature)
  for (i in 1:length(idxERCC)){
    dat$Ratio[i] <- as.character(expDat$idColsAdj$Ratio)[match(dat$Feature[i],
                                           expDat$idColsAdj$Feature)]
  }

  #idx1 <- match( expDat$idColsAdj$Feature, dat$Feature,
  #               nomatch=F)
  #dat$Ratio[idx1] <- as.character(expDat$idColsAdj[idx1,c(4)])
  dat$Ratio <- as.factor(dat$Ratio)
  
  
  
  myYLim = plotInfo$myYLimMA
  myXLimMA = plotInfo$myXLimMA
  
  filenameRoot = sampleInfo$filenameRoot
  sample1 = expDat$sampleNames[1]
  sample2 = expDat$sampleNames[2]
 
  
  #idCols = expDat$idColsAdj
  r_m.res = expDat$Results$r_m.res
  
  #cutoffs = LODR.annot.ERCC$cutoffs
  cutoffs = expDat$Results$LODR.annot.ERCC$countCutoffs
  
  FCcode = erccInfo$FCcode
  legendLabels=erccInfo$legendLabels
  
  
  spikeFraction = expDat$spikeFraction
  
  r_m.mn = exp(r_m.res$r_m.mn)
  r_m.UL = exp(r_m.res$r_m.upper)
  r_m.LL = exp(r_m.res$r_m.lower)
  
  #theme_update(legend.justification=c(1,0), legend.position=c(1,0))
  
  
  colScale <- plotInfo$colScale
  fillScale <- plotInfo$fillScale
  
  
  if (length(cutoffs) > 0){
    cat("\nLODR estimates are available to code ratio-abundance plot\n")  
  }
  
  maStats <- function(x, c1, c2){
    c(mean(log2(x[c1])-log2(x[c2])),sd(log2(x[c1])-log2(x[c2])),log2(mean(x)))
  } 
  
  totCol <- ncol(dat[-c(1:2)])
  
  if(odd(totCol)) stop("Uneven number of replicates for the two sample types")
  
  maStatDat <- data.frame(t(apply(dat[-c(1:2)],1,maStats, 
                                  c1 = c(1:(totCol/2)),
                                  c2 = c(((totCol/2)+1):totCol))))
  colnames(maStatDat) <- c("M.Ave","M.SD","A")
  maDatAll <- cbind(dat, maStatDat)
  
  maDatAll <- maDatAll[which(is.finite(maDatAll$M.Ave)),]

  ### Now subset and continue with just ERCCs
  maData <- subset(maDatAll, Ratio != "Endo")
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
  }else{
    maData$Empirical = maData$Nominal
  }
  
  
  
  cutERCCs = unique(maData$Feature[which(is.na(maData$M.SD))])
    
  if(length(cutERCCs) != 0){
    cat(paste("\nThese ERCCs were not included in the ratio-abundance plot, \n",
              "because not enough non-zero replicate measurements of these \n",
              "controls were obtained for both samples:\n"))
    
    cat(paste(as.character(cutERCCs),collapse="\n"))
    
    maData = maData[-(which(maData$Feature %in% cutERCCs)),]
    searchCut = paste(cutERCCs,collapse="|")
    maData = maData[-(grep(searchCut,maData$Feature)),]
  }
    
    #maDataAveSD$Feature = as.factor(as.character(maDataAveSD$Feature))  
    maData$Feature = as.factor(as.character(maData$Feature))
    #write.csv(maData,file = paste(filenameRoot,"maDataFinite.csv",sep = "."))
  
  #avexlabel = expDat$ERCCxlabelAve
  if(sampleInfo$datType == "count"){
    avexlabel = "Log2 Average of Read Depth Normalized Counts"
    ymalabel = "Log2 Ratio of Read Depth Normalized Counts"
  }
  if(sampleInfo$datType == "array"){
    avexlabel = "Log2 Average of Normalized Intensity"
    ymalabel = "Log2 Ratio of Normalized Intensity"
   # myXLimMA = c(min(maData$A)-1, max(maData$A)+1)
   # expDat$plotInfo$myXLimMA <- myXLimMA
  }
  xlabel = xlab(avexlabel)
  
    # Estimate SD for all ratio measurements
    # For the ERCCs that are plotted take the log ratio data and subtract the log 
    # nominal ratio, this will normalize the data
    normRats = maData$M.Ave - log2(maData$Nominal)
    #print(normRats)
    
    # Take the sd of all of the normalized log ratio data to find a global SD for the ERCCs at this site
    sdGlobal = sd(normRats)
    #cat("\nGlobal Ratio SD for this sample pair is: ")
    #cat(sdGlobal, "\n")
    
    ratVarDat <- maData
    #ratVarDat <- subset(ratVarDat, A > 0)
    #ratVarDat$A <- log2((2^(ratVarDat$A.Ave))/(spikeFraction))
    
    sdRatioplot <- ggplot(ratVarDat) + geom_point(data = ratVarDat, aes(x = A,
                                                               y = M.SD, 
                                                               colour = Ratio),
                                         size = 5, alpha = alphaPoint)+colScale 
      
    
    lmfit <- lm(formula = M.SD ~ 1, data = ratVarDat, weights = 1/((M.SD)^2))
    
    stdevCoef <- signif(summary(lmfit)$coefficients, digits = 4)
    
    params<-list(Mo = 0.1, lambda = 2) # setup default params.
    f = coef(lmfit)[[1]]  
    #print("Using default parameters")
    #print(params)
    #print(f)
    
    while(TRUE){
      
      expmod<-NULL
      # does not stop in the case of error
      try(expmod<-nls(formula=M.SD ~ Mo*exp(-lambda * A) + f,
                      data = ratVarDat,
                      start = params), silent = TRUE)
      if(!is.null(expmod))break; # if nls works, then quit from the loop
      #print("Try changing initial guess for parameters")
      params<-list(Mo = sdGlobal, lambda = 0.2)# change the params for nls
      #print(params)
      
    }
    
    if(!is.null(expmod)){
      stdevCoef <- rbind(stdevCoef, signif(summary(expmod)$coefficients, digits = 4))
    } 
    
    #print(stdevCoef)
    row.names(stdevCoef)<- c(expression("V[min]"), 
                             expression("V[max]"),
                             "lambda")
    #print(stdevCoef[,c(1:2)])
    stdevTable <- stdevCoef[c(1:3),c(1:2)]
    colnames(stdevTable)<- c(expression("Estimate"), expression("`Standard Error`"))
    
    # add fitted curve
    lineDat <- data.frame("xdat" = ratVarDat$A,"ydat" = 
                            predict(lmfit, x = ratVarDat$A)) 
    if(!is.null(expmod)){
#       fitDat <- data.frame("xdat" = ratVarDat$A,
#                            "ydat" = predict(expmod,
#                                             list(x = ratVarDat$A))) 
      xdat = seq(myXLimMA[1],myXLimMA[2],length.out=50)
      fitDat <- data.frame("xdat" = xdat,
                           "ydat" = predict(expmod,newdata=
                                            list(A = xdat))) 
    } 
  
  sdRatioplotFit <- sdRatioplot +# geom_line(data = lineDat,
                                  #            aes(x=xdat,y=ydat)) +
    ylab(paste0("sd(", ymalabel, ")")) + xlabel +
    coord_cartesian(xlim = myXLimMA, ylim = c(0,1.5*max(ratVarDat$M.SD))) +
    annotation_custom(grob=tableGrob(stdevTable,parse=T, gpar.corefill = 
                                  gpar(fill = "grey85", col = "white"),
                                gpar.rowfill = gpar(fill = "grey80",
                                                    col = "white"),
                                gpar.colfill = gpar(fill = "grey80",
                                                     col = "white")),
                      xmin = 0.25*max(ratVarDat$A),xmax = max(ratVarDat$A),
                      ymax = max(ratVarDat$M.SD),
                      ymin = 0.75*max(ratVarDat$M.SD)) +
    #annotate("text", x = 0.5*max(ratVarDat$A), y = 0.5*max(ratVarDat$M.SD),
    #         label = paste0("y = V[max]","e^{-lambdax}",
    #                                        " + V[min]")) +
    theme_bw() + theme(legend.justification=c(0,0), legend.position=c(0,0)) 
    
    
    if(!is.null(expmod)){
      sdRatioplotFit <- sdRatioplotFit + 
        geom_line(data = fitDat, aes(x = xdat, y = ydat)) + xlabel      
    } 
    
    if (length(cutoffs)>0){
      maData$LODR = "below"
      FCcodeC = FCcode[-c(which(FCcode$FC == 1)),]
      for (i in 1:length(cutoffs)){
        maData$LODR[which((maData$A > cutoffs[i])&
                                 (maData$Ratio ==
                                    FCcodeC$Ratio[i]))] = "above"
      }
      
      maData$LODR <- as.factor(maData$LODR)
      rm_dat = data.frame("Mean r_m" = signif(r_m.mn,digits = 4),
                          "Lower Limit" = signif(r_m.LL,digits = 4), 
                          "Upper Limit" = signif(r_m.UL,digits = 4))
      colnames(rm_dat) <- c("Mean", expression(paste("95% CI Lower Bound")), 
                            expression(paste("95% CI Upper Bound")))
      rownames(rm_dat) <- expression(r[m])
      
      # Plot ratio-signal data coding points below the LODR with open circles
      maPlot <- ggplot(maData, aes(x = A, y = M.Ave) ) + 
        geom_point(data = subset(maDatAll, maDatAll$Ratio == "Endo"), aes(x = A,
                                                                  y = M.Ave),
                   colour = "grey80", alpha = 0.5) +
        geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD, 
                          colour = Ratio),size = 1,alpha = alphaPoint) + 
        geom_point(aes(colour = Ratio),size = 5, alpha = alphaPoint) +
        geom_point(data = subset(maData, (LODR == "below")),
                   colour = "white",size = 2.5) + 
        geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                   linetype = "longdash",alpha = 0.7) + 
        ylab(ymalabel) + xlabel + 
        coord_cartesian(xlim = myXLimMA, ylim = myYLim) + colScale + 
        annotation_custom(tableGrob(rm_dat,parse=T, 
                                    gpar.corefill = gpar(fill = "grey85",
                                                        col = "white"), 
                                    gpar.rowfill = gpar(fill = "grey80",
                                                        col = "white"),
                                    gpar.colfill = gpar(fill = "grey80",
                                                        col = "white")), 
                          #xmin = quantile(maData$A,probs=0.25),
                          #xmax = max(maData$A),
                          ymin = (myYLim[2]) - 0.25*myYLim[2], 
                          ymax = myYLim[2]) + 
        scale_y_continuous(breaks = seq(myYLim[1],myYLim[2],1))+ theme_bw()+
        theme( legend.justification = c(1,0),legend.position=c(1,0))
        
    }else{
      # Plot ratio-signal data without LODR coding
      maPlot <- ggplot(maData, aes(x = A, y = M.Ave, colour = Ratio) ) + 
        geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD), 
                      alpha = 0.3) + geom_point(size = 5, alpha = alphaPoint) + 
        geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                   linetype = "longdash",alpha = 0.7) + 
        ylab(ymalabel) + xlabel + colScale+ 
        coord_cartesian(xlim = myXLimMA, ylim = myYLim) + theme_bw()+
        theme( legend.justification = c(1,0),legend.position=c(1,0))
    }
  
  ratioVarPlot <- ggplot(maData) + geom_violin(aes(x = Ratio, 
                                                   y = M.SD, 
                                                   fill = Ratio), 
                                               alpha = alphaPoint) +
  ylab("SD of Log2 Ratios") + colScale + fillScale + theme_bw()

  
  expDat$Figures$maPlot <- maPlot
  #expDat$Results$ratVarDat <- ratVarDat
  expDat$Results$modRatVar <- stdevCoef
  expDat$Results$maDatAll <- maDatAll
  #expDat$Figures$ratioSDPlot <- sdRatioplotFit
  expDat$Figures$ratioSDPlot <- ratioVarPlot
  
  return(expDat)
}