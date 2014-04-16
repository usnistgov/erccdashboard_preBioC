#' Estimate Limit of Detection of Ratios (LODR)
#'
#' @param exDat    list, contains input data and stores analysis results
#' @param kind      "ERCC" or "Sim"
#' @param prob      probability, ranging from 0 - 1, default is 0.9
#' 
#' @details
#' This is the function to estimate a limit of detection of ratios (LODR) for a
#' a chosen probability and threshold p-value for the fold changes in the ERCC
#' control ratio mixtures.
#' @export
#' 
estLODR <- function(exDat,kind = "ERCC", prob=0.9){
  ##############
  #### LODR ####
  ##############
  ### Select fitting parameters
  # pval.cutoff default is .001
  # probability cutoff is .9
  # kind can be ERCC or Sim
  sampleInfo <- exDat$sampleInfo
  erccInfo <- exDat$erccInfo
  plotInfo <- exDat$plotInfo
  
  pval.cutoff<- exDat$Results$p.thresh 
  filenameRoot<-sampleInfo$filenameRoot
  FCcode <- erccInfo$FCcode
  #FCcode = data.frame(Ratio = c("4:1","1:1","1.5:1","2:1"),FC =  c(4,1,.667,.5))
  legendLabels = erccInfo$legendLabels
  #legendLabels = FCcode$Ratio

  colScale <- plotInfo$colScale
  fillScale <- plotInfo$fillScale
  
  ### Read in results for LODR ####
  ### Modify this file path to work with other ERCC data
  ### Note that you must format your provided data in the same manner as the provided example
  ### i.e. the column headings of your selected datamust match those of the provided examples
  if(file.exists(paste(filenameRoot,kind,"Pvals.csv"))){
    pval.res<-read.csv(paste(filenameRoot,kind,"Pvals.csv"),header=TRUE)  
  }else{
    cat(paste0("\n",filenameRoot," ",kind," Pvals.csv is missing, ",
              "because simulated DE data was not generated \nand tested for ",
              "differential expression. Continuing with analysis...\n"))
    return(exDat)
  }
   
  names(pval.res)[1]= "Feature"
  
  #library(locfit)
  #library(gridExtra)
  cat(paste("\nEstimating",kind,"LODR\n"))
  #### Function used to estimate LODR and its uncertainty
  LODR<-function(pval,mn,cutoff,prob){
    cutoff<-log10(cutoff)
    fit<-locfit(log10(pval)~lp(log10(mn)))
    X<-preplot(fit,band="pred",newdata=log10(mn))
    
    ### See plot of fit
    #plot(fit,band="pred",get.data=TRUE,xlim=range(log10(mn)))
    
    find.mn<-function(mn,fit,cutoff,prob){
      X<-preplot(fit,newdata=mn,band="pred")
      (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
    }
   
    rng.mn<-range(log10(mn))
    
  #  ### Search in sections to get first crossing 
   # segmented.search<-function(fit){
# 	ppp<-.2
#     t.lodr<-optimize(find.mn,c(rng.mn[1],sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
#     while(t.lodr$objective>.0001&ppp<=1){
#         t.lodr<-optimize(find.mn,c(sum(rng.mn*c(1-ppp+.2,ppp-.2)),sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
#         ppp<-ppp+.2
#     }
#     t.lodr
# }   
  ### Search in sections to get first crossing 
  segmented.search<-function(fit){
    X<-preplot(fit,newdata=min(log10(mn)),band="pred")
    if((X$fit+qnorm(prob)*X$se.fit)<cutoff){
      t.lodr<-list(minimum=min(log10(mn)),objective=0)
    } else{
      ppp<-.2
      t.lodr<-optimize(find.mn,c(rng.mn[1],sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
      while(t.lodr$objective>.0001&ppp<=1){
        t.lodr<-optimize(find.mn,c(sum(rng.mn*c(1-ppp+.2,ppp-.2)),sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
        ppp<-ppp+.2
      }
    }
    t.lodr
  }    
  
  lodr<-segmented.search(fit)
     
    ### Bootstrap to estimate uncertainty
    lodr.boot<-NULL
    for(ii in 1:500){
      y.boot<-X$fit+sample(residuals(fit),length(mn))
      #if(ii %in%c(20*(1:5))){points(log10(mn),y.boot,col=j); j<-j+1}
      fit.boot<-locfit(y.boot~lp(log10(mn)))
      lodr.boot<-c(lodr.boot,segmented.search(fit.boot)$minimum)
      if (ii %% 100 == 0){
        cat ("...")
      }
    }
    return(c(lodr$objective,lodr$minimum,quantile(lodr.boot,c(.05,.95))))
  }
    

  #### Apply to loaded data ####
  lodr.resPlot<-NULL; set.seed(1)
  lodr.resLess<-NULL; set.seed(1)
  
  #lodr.res<-NULL; set.seed(1)
  lineDat <-NULL
  
  pval.res$Ratio = "Ratio"
  for (i in 1:nlevels(FCcode$Ratio)){
    pval.res$Ratio[which(pval.res$Fold == FCcode$FC[i])] = as.character(FCcode$Ratio[i])
  }
  
  #pval.res$Ratio <- as.factor(pval.res$Ratio)
  pval.res$Ratio <- factor(pval.res$Ratio, levels = as.character(FCcode$Ratio))
  for(i in 1:nlevels(FCcode$Ratio)){  ### Analyze each fold change subgroup separately
    
    grp<-which(pval.res$Ratio==FCcode$Ratio[i])  ### identify subgroup members

    x<-pval.res$MnCnt[grp][pval.res$Pval[grp]!=0]; y<-pval.res$Pval[grp][pval.res$Pval[grp]!=0]
    
    fit<-locfit(log10(y)~lp(log10(x)))
    #x.new<-log10(pval.res$MnCnt[grp])
    x.new<-seq(min(log10(x)),max(log10(x)),length.out=100)
    X<-preplot(fit,band="pred",newdata=x.new)

    x.new<-10^x.new
    fitLine = 10^(X$fit)
    fitUpper = 10^(X$fit+qnorm(prob)*X$se.fit)
    fitLower = 10^(X$fit-qnorm(prob)*X$se.fit)
    fitData = data.frame(x.new,fitLine,fitUpper,fitLower,Ratio = FCcode$Ratio[i])
    lineDat = rbind(lineDat,fitData)

# ### Estimate LODR
#lodr.res<-rbind(lodr.res,c(abs(log2(FCcode$FC[i])),LODR(y,x,cutoff=pval.cutoff,prob=prob)))

    ### Estimate LODR
    if(FCcode$FC[i]!=1){
      t.res<-LODR(y,x,cutoff=pval.cutoff,prob=prob)
      t.res[-1]<-signif(10^t.res[-1],2)
      if(t.res[1]>.01){
        t.res[2]<-Inf
        t.res[3:4]<-NA
      }
      t.resLess <- t.res
      t.resLess[-1][t.resLess[-1]==signif(min(x),2)]<-paste("<",signif(min(x),2),sep="")
      t.res[-1][t.res[-1]==signif(min(x),2)] <- Inf
    }
    if(FCcode$FC[i]==1){
      t.res<-rep(NA,4)
      t.resLess<-rep(NA,4)
    } 
    lodr.resPlot<-rbind(lodr.resPlot,c(round(abs(log2(FCcode$FC[i])),3),t.res))
    lodr.resLess <- rbind(lodr.resLess,c(round(abs(log2(FCcode$FC[i])),3),t.resLess))
    
}
  pval.res = subset(pval.res,pval.res$Pval != 0)
  #print(lodr.resLess)
  #print(lodr.resPlot)
  
  colnames(lodr.resLess)[1:3]<-c("|log2(Fold)|","MinError","Estimate")
  #print(paste(filenameRoot,kind,"LODR Estimates for P-value Threshold:",pval.cutoff,"; Probability:",prob))
  
  colnames(lodr.resPlot)[1:3]<-c("Ratio","MinError","Estimate")
  colnames(lodr.resLess)[1:3]<-c("Ratio","MinError","Estimate")
  lodr.resPlot <- as.data.frame(lodr.resPlot)
  lodr.resLess <- as.data.frame(lodr.resLess)
  lodr.resPlot$Ratio <- as.character(legendLabels)
  lodr.resLess$Ratio <- as.character(legendLabels) 
  annoTable = lodr.resLess[-c(2)]
  colnames(annoTable) <- c("Ratio",expression("LODR Estimate"), expression("90% CI Lower Bound"), expression("90% CI Upper Bound"))
  annoTable <- annoTable[-which(annoTable$Ratio == "1:1"),]
  cat("\n")
  print(annoTable)
  
  arrowDat = data.frame(Ratio = FCcode$Ratio, FC = FCcode$FC, x = lodr.resPlot[,3], y = pval.cutoff, xend = lodr.resPlot[,3], yend = 0)
  arrowDat$x[grep('<',lodr.resLess[,3])] <- Inf
  arrowDat = arrowDat[-which(arrowDat$FC == "1"),]
  arrowDat = arrowDat[which(is.finite(arrowDat$x)),]
  
  if(dim(arrowDat)[1] == 0){
    cat(paste("\nError! Estimated distribution of p-values does not cross \n",
              "threshold p-value, may be due to insufficient data quantity\n"))
    break
  }
  if(exDat$sampleInfo$datType == "array"){
    xlabDE <- xlab("Average Fluorescence Intensity")
    xrange <- c(min(pval.res$MnCnt),max(pval.res$MnCnt))
    legendPos <- theme(legend.justification=c(1,0), legend.position=c(1,0))
  }
  if(exDat$sampleInfo$datType == "count"){
    xlabDE <- xlab("Average Counts")
    xrange <- c(1,max(pval.res$MnCnt))
    legendPos <- theme(legend.justification=c(0,0), legend.position=c(0,0))
  }
  
  ## create inset table
  my_table <- tableGrob(d=annoTable,
                        show.rownames=F,
                        gpar.coretext =gpar(fontsize=14),
                        gpar.coltext=gpar(fontsize=14),
                        gpar.rowtext=gpar(fontsize=14),
                        gpar.corefill = gpar(fill = "grey85", col = "white"), 
                        gpar.rowfill = gpar(fill = "grey80", col = "white"),
                        gpar.colfill = gpar(fill = "grey80", col = "white"))
  
#print(min(log(pval.res$Pval)))
#   if (((min(pval.res$LogPval)) < log(1e-15)) & (yrange[1]>=1e-15)){
#     cat("\nVery small p-values, adjusting y-limits to show all data\n")
#     ymin <- min(pval.res$Pval)/1e5
#     yrange <- c(ymin, 1e1)
#   }
  LODRplot <- ggplot(pval.res, aes(x = MnCnt, y = Pval, colour = Ratio)) + 
    geom_point(size = 6) + 
    scale_x_log10(limits = xrange) + 
    #scale_y_log10(breaks = c(1e-12,1e-10,1e-8,1e-5,1e-4,1e-3,1e-2,1e-1,1e0), 
    #              limits = yrange) + 
    scale_y_log10(breaks = c(1e-12,1e-10,1e-8,1e-5,1e-4,1e-3,1e-2,1e-1,1e0))+
    geom_ribbon(data = lineDat, aes(x = x.new, 
                                    y = fitLine, 
                                    ymin=fitLower, 
                                    ymax=fitUpper,
                                    fill = Ratio), alpha = 0.3,colour = NA,show_guide =F) + 
    geom_line(data = lineDat,aes(x = x.new,
                                 y=fitLine, 
                                 colour = Ratio),show_guide = F) + 
    colScale + fillScale + xlabDE + ylab("DE Test P-values") + 
    geom_hline(yintercept = pval.cutoff, linetype = 2, size = 2 ) + 
    geom_segment(data = arrowDat, aes(x = x,
                                      y = y, 
                                      xend = xend,
                                      yend = yend, 
                                      colour = Ratio), 
                 lineend = "round", arrow = arrow(length =unit(0.5,"cm")), 
                 size = 2, alpha = 0.6) + theme_bw()+ legendPos
    
    
  ### final result 
  
#   Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2, 0.5), c("null", "null")))
#   #grid.show.layout(Layout)
#   vplayout <- function(...) {
#     grid.newpage()
#     pushViewport(viewport(layout = Layout))
#   }
#   
#   subplot <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#   
#   mmplot <- function(a, b) {
#     vplayout()
#     print(a, vp = subplot(1, 1))
#     print(b, vp = subplot(2, 1))
#   }
#   
# 
#   annotLODRplot <- mmplot(LODRplot, arrangeGrob(my_table))
  #vplayout()
  
  annotLODRplot <- arrangeGrob(LODRplot, arrangeGrob(my_table), 
                               ncol = 1, heights = c(2,0.5))
  #annotLODRplot <- grid.arrange(LODRplot, arrangeGrob(my_table), ncol = 1, heights = c(2,0.5))
  #print(annotLODRplot)
  
  nam <- paste0("lodr",kind,"Plot")
  exDat$Figures$plotLODR <- annotLODRplot
  #print(exDat$Figures$plotLODR)
  names(exDat$Figures)[which(names(exDat$Figures) == "plotLODR")] <- nam
  
  nam <- paste("lodr.res",kind,sep = ".")
  exDat$Results$lodr.res <- lodr.resLess
  names(exDat$Results)[which(names(exDat$Results) == "lodr.res")] <- nam
  
  return(exDat)
}