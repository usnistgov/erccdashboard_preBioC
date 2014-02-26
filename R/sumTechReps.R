sumTechReps <- function(expDat,  libeList = libeList){
  if(is.data.frame(expDat$expressDatSumNorm)){
    print ("Tech Reps already summed")
    break
  }
  sampleInfo = expDat$sampleInfo
  expressDat = expDat$TranscriptsAB 
  designMat = expDat$designMatAB 
  sampleNameList = expDat$sampleNames 
  
  platform <- sampleInfo$platform
  siteName <- sampleInfo$siteName
  libeSizeNorm <- sampleInfo$libeSizeNorm
  totalSeqReads <- sampleInfo$totalSeqReads
  totalReads <- expDat$totalReads
  
  myDataERCC = expressDat[-c(1)]
  expressDatSum = expressDat[c(1)]
  totalReadSum = NULL
  
  #if(type == "sum"){
    if(libeSizeNorm == T){
      if (totalSeqReads == T){
        totalReads = cbind(totalReads,designMat)
        #print(head(totalReads))
        #totalReadSum = NULL
      }
      # In a nested for loop sum the A-D lanes and flow cells for each library
      for (sampleName in 1:length(sampleNameList)){
        for (libeNum in 1:length(libeList)){
          select = subset(designMat, (Sample == sampleNameList[sampleName])&(Library == libeList[libeNum]))
          select <- as.data.frame(lapply(select,as.character))
          select <- as.data.frame(lapply(select,as.factor))
          
          expressDatForSum = myDataERCC[c(match(select$countSet, names(myDataERCC)))]
          sumCounts = rowSums(expressDatForSum)
          colName = paste("SEQC",platform,siteName,sampleNameList[sampleName],libeNum, sep = "_")
          expressDatSum = cbind(expressDatSum,sumCounts)
          names(expressDatSum)[ncol(expressDatSum)] = colName
          
          if(totalSeqReads == T){           
            totalReadsSummed = totalReads$totalReads[c(match(select$countSet,totalReads$countSet))] 
            sumReads = sum(totalReadsSummed)
            totalReadSum = c(totalReadSum, sumReads)
          }   
        }
      }
#### Commented 131101 don't need to do filtering here, needs to be in 
      # loadExpMeas file for users
#       dat = expressDatSum[-c(1)]
#       ## Filter out genes with average counts less than 1 or fewer than three libraries with at least 1 count
#       expressDatSum<-subset(expressDatSum, (rowMeans(dat)>1)&(rowSums(dat!=0)>2)); 
#       
      # create Data frame of the Sample library pairs to return to the workspace for interlab analysis
      allSampleLibeData = expressDatSum 
      
      datCols = expressDatSum[-c(1)]
      
      if (totalSeqReads == F){
        TranscriptsAll = expressDatSum  
        TranscriptMappedReadSums = colSums(TranscriptsAll[-c(1)],na.rm = T)
        libeSize = TranscriptMappedReadSums
        libeSize = libeSize/(10^6) #per million mapped reads
        print(libeSize)
      }else{
        TranscriptsAll = expressDat[-c(grep("ERCC-0", expressDat$Feature)),]  
        libeSize = totalReadSum
        libeSize = libeSize/(10^6) #per million sequenced reads
      }
  
      #Library size normalize the data  
      libAdjust = sweep(datCols, 2, libeSize,"/")
      expressDatSummed = cbind(expressDatSum[c(1)], libAdjust)
    }
    else{
      # In a nested for loop sum the A-D lanes and flow cells for each library
      for (sampleName in 1:length(sampleNameList)){
        for (libeNum in 1:length(libeList)){
          select = subset(designMat, (Sample == sampleNameList[sampleName])&(Library == libeList[libeNum]))
          select <- as.data.frame(lapply(select,as.character))
          select <- as.data.frame(lapply(select,as.factor))
          
          expressDatForSum = myDataERCC[c(match(select$countSet, names(myDataERCC)))]
          meanCounts = rowMeans(expressDatForSum)
          #meanCounts = rowSums(expressDatForSum)
          colName = paste("SEQC",platform,siteName,sampleNameList[sampleName],libeNum, sep = "_")
          expressDatSum = cbind(expressDatSum,meanCounts)
          names(expressDatSum)[ncol(expressDatSum)] = colName
        }
      }
      dat = expressDatSum[-c(1)]
      ## Filter out genes with average counts less than 1 or fewer than three libraries with at least 1 count
      expressDatSum<-subset(expressDatSum, (rowMeans(dat)>1)&(rowSums(dat!=0)>2)); 
      
      # create Data frame of the Sample library pairs to return to the workspace for interlab analysis
      allSampleLibeData = expressDatSum 
      expressDatSummed = expressDatSum
    }
    
    expDat$expressDatSumNorm = expressDatSummed
    expDat$expressDatSumNoNorm = allSampleLibeData
    expDat$totalReadSum = totalReadSum
    if (totalSeqReads == T) expDat$totalReads = totalReadSum
    expDat$designMatSum = getDesignMat(expDat$expressDatSumNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')
    return(expDat)
    
  #}
#   if(type == "mean"){
#     myDataERCC = expressDat[-c(1)]
#     expressDatMean = expressDat[c(1)]
#     
#     # In a nested for loop sum the A-D lanes and flow cells for each library
#     for (sampleName in 1:length(sampleNameList)){
#         select = subset(designMat, (Sample == sampleNameList[sampleName]))
#         select <- as.data.frame(lapply(select,as.character))
#         select <- as.data.frame(lapply(select,as.factor))
#         expressDatForMean = myDataERCC[c(match(select$countSet, names(myDataERCC)))]
#         meanCounts = rowMeans(expressDatForMean)
#         colName = paste("SEQC",platform,siteName,sampleNameList[sampleName], sep = "_")
#         expressDatMean = cbind(expressDatMean,meanCounts)
#         #print(head(expressDatMean))
#         names(expressDatMean)[ncol(expressDatMean)] = colName
#     }
#     
#     # create Data frame of the Sample library pairs to return to the workspace for interlab analysis
#     expDat$expressDatCombinedMean = expressDatMean
#     #expDat$allSampleLibeData = NULL
#     #expDat$totalReadSum = NULL
#     
#     return(expDat)
#   }
}