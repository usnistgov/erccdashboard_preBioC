#loadArray
### get the COH array
  fileName ="/Users/smunro/Documents/NISTMunro/Projects/erccdashboard201309/data/QCmicroarrays_COH_UTSW/COHdata110909.txt"){
  COHarray <- read.table(fileName,header = T)
  COHarrayAB <- COHarray[-c(11:16)]
  y <- COHarrayAB[c(1,2,5:10)]
  y$PROBE_ID<- as.character(y$PROBE_ID)
  idxERCC <- grep("ERCC-",y$TargetID)
  print(head(y))
  # copy the ERCC target ids over to the probe ID column
  y$PROBE_ID[idxERCC] <- as.character(y$TargetID[idxERCC])
  y$PROBE_ID <- as.factor(y$PROBE_ID)
  y$TargetID<- NULL
  colnames(y) <- c("Feature", "UHRR_3","UHRR_2","UHRR_1","HBRR_3","HBRR_2","HBRR_1")

  
  ### get the UTSW array
  fileName ="/Users/smunro/Documents/NISTMunro/Projects/erccdashboard201309/data/QCmicroarrays_COH_UTSW/UTSWdata110909.txt"
    array <- read.table(fileName,header = T)
    arrayAB <- array[c(1,2,4,8,12,5,9,13)]
    
    y <- arrayAB
  
    array1 <- as.character(y$TargetID)
    array2 <- as.character(y$ProbeID)
  
    Feature <- paste(array1,array2,sep="_")
    
    idxERCC <- grep("ERCC-",y$TargetID)
    Feature[idxERCC] <- as.character(y$TargetID[idxERCC])
  
    y$TargetID <- Feature
    y[c(2)]<- NULL
  
    
    print(head(y))
    colnames(y) <- c("Feature", "UHRR_1","UHRR_2","UHRR_3","HBRR_1","HBRR_2","HBRR_3")
    