compareReplicates <- function(repDat = expressDat, idCols, siteName, select, repType = "Lane", SampleName = "A"){
  theme_set(theme_bw(base_size=16))
  theme_update(legend.justification=c(1,0), legend.position=c(1,0))
  select <- as.data.frame(lapply(select,as.character))
  select <- as.data.frame(lapply(select,as.factor))
  repDatDist = repDat[c(1,match(select$countSet, names(repDat)))]
  repDatDistERCCs = merge(idCols[c(1,4)],repDatDist)
  meanDist = apply(X = repDatDistERCCs[-c(1:2)],MARGIN=1,FUN=mean)
  sdDist = apply(X = repDatDistERCCs[-c(1:2)],MARGIN=1,FUN=sd)
  
  expressDisttoMean = repDatDistERCCs[,c(1,2)]
  lastCol = ncol(repDatDistERCCs)
  for (i in c(3:lastCol)){
    DisttoMean = data.frame(repDatDistERCCs[c(i)]/meanDist)
    colName = paste( repType,as.character(i-2), sep = "_")
    names(DisttoMean)[1] = colName
    expressDisttoMean = cbind(expressDisttoMean,DisttoMean )  
  }
  head(expressDisttoMean)
  
  expressDisttoMean_l = melt(expressDisttoMean, id.vars=c("Feature", "Ratio") )
  names(expressDisttoMean_l)[3]= "Replicate"
  
  # plot distributions as function of library
  pairedDensity = ggplot(expressDisttoMean_l) + geom_density(aes(x = value, fill = Replicate),alpha = 0.35,linetype = 0) + coord_flip() + geom_rug(aes(x = value, colour = Replicate))  + xlab(paste("L Distributions Per",repType)) + xlim(-2,2)
  
  meanDist = data.frame("L.value" = apply(X = expressDisttoMean[c(3:ncol(expressDisttoMean))],MARGIN=2,FUN=mean,na.rm = T))
  meanDist$Dist = row.names(meanDist)
  row.names(meanDist) <- NULL
  
  kruskal.res = kruskal.test(as.vector(expressDisttoMean[c(3:ncol(expressDisttoMean))]))
  print("Kruskal-Wallis Test P-value")
  print(kruskal.res$p.value)
  if(kruskal.res$p.value < 0.01){
    print("Kruskal-Wallis Test P-value is Less Than 0.01, Use Wilcoxon Test to Compare Distributions")
    #as.vector(expressDisttoMean[c(3:ncol(expressDisttoMean))]
    for (i in 3:ncol(expressDisttoMean)){
      print(names(expressDisttoMean)[i])
      singleDist = expressDisttoMean[,i]
      groupDist = as.matrix(expressDisttoMean[,-c(1,2,i)])
      groupDist_l = as.vector(groupDist)
      wilcox.res = wilcox.test(x = singleDist,y = groupDist)
      print(wilcox.res$p.value)
      # build a rank ordered data frame with a replicate column and a p-value column
      
    }
  }
  else{
    print("Cannot reject the null hypothesis that distributions are similar")
  }
  
  return(list(pairedDensity,meanDist))
}