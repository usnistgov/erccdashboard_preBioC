loadExpMeasFromRaw <- function(expDat){
  sampleInfo = expDat$sampleInfo
  # Pull idCols out of expDat
  idCols <- sampleInfo$idColsSRM 
  
  # Import data based on analysis type from SEQC main project
 
 if (sampleInfo$study == "SEQC_Main"){
   nameKey <- data.frame(oldName = c("AGR","BGI","CNL","COH","MAY","NVS","NWU",
                                     "PSU","SQW"), newName = c("Lab1","Lab2",
                                                               "Lab3", "Lab4",
                                                               "Lab5", "Lab6",
                                                               "Lab7", "Lab8", 
                                                               "Lab9"))
   oldName <- nameKey$oldName[which(nameKey$newName == sampleInfo$siteName)]
   sampleInfo$sample1Name <- "A"
   sampleInfo$sample2Name <- "B"
   
   # Get total Reads (from the original fastq files)
   if(sampleInfo$totalSeqReads == T){
     if(sampleInfo$platform != "ROC"){
       mainReads = read.delim(paste("Data/MainSEQC/all_qc_results/",
                                    sampleInfo$platform,"_",oldName,"_qc_results.txt",
                                    sep=""))
       mainReads = mainReads[c(1:9)]
       if(sampleInfo$analysis == "NIST"){
         if (oldName == "NWU") mainReads$sampleName <- 
           sub("_NWU00001.*","", mainReads$sampleName)
         if (oldName == "SQW") mainReads$sampleName <- 
           sub("_23311016_20120217.*","", mainReads$sampleName)
         if (oldName == "PSU") mainReads$sampleName <- 
           sub("_23302023_20120308.*","", mainReads$sampleName) 
       }  
     }else{
       mainReads = read.delim("Data/MainSEQC/RocheResults/Mappability.txt",
                              sep = " ")
       mainReads = mainReads[c(1,2,4,5)]
       colnames(mainReads) <- c("Site","Sample", "Replicate","TotalReads") 
     }
   }
   if (sampleInfo$analysis == "NCTR"){
     dataFile = paste("Data/MainSEQC/SEQC_MAIN_ILM_rawCounts_ZSu.2012_04_14/",
                      "SEQC_MAIN_ILM_",oldName,"_TranscriptCounts_ZSU.txt",
                      sep = "")
     Transcripts = read.delim(dataFile)
   }
   if (sampleInfo$analysis == "NIST"){
     dataFile = paste("Data/MainSEQC/SEQC_LifescopeCounts/SEQC_LIF_",
                      oldName,"_LifeScope.csv", sep = "")
     Transcripts = read.csv(dataFile, header = T)
     newColnames <- colnames(Transcripts)
     #colnames(Transcripts) <- sub("_F3.*","", newColnames)
     if (oldName == "NWU") colnames(Transcripts) <- 
       sub("_NWU00001.*","",newColnames)
     if (oldName == "SQW") colnames(Transcripts) <- 
       sub("_23311016_20120217.*","", newColnames)
     if (oldName == "PSU") colnames(Transcripts) <- 
       sub("_23302023_20120308.*","", newColnames)
   } 
   if(sampleInfo$analysis == "WEHI"){
        dataFile = "Data/MainSEQC/RocheResults/RefSeq-All-Genes.txt"
        Transcripts = read.delim(dataFile, header = T)
        dataAB = Transcripts[-c(1)]
        Transcripts =cbind(Transcripts[c(1)], dataAB[order(colnames(dataAB))])
        print(head(Transcripts))
   }
   
   # force names to be ERCC- and first column name to Feature
   names(Transcripts)[1] = "Feature"
   Transcripts$Feature = gsub("ERCC_","ERCC-",Transcripts$Feature)
   Transcripts$Feature = gsub(":Gene_AceView08","",Transcripts$Feature)
   Transcripts$Feature = gsub(":Gene_RefSeq","",Transcripts$Feature)

   # get data frames with just the ERCCs and just the human genes
   TranscriptsERCCOnly = Transcripts[c(grep("ERCC-0", Transcripts$Feature)),]
   TranscriptsHumanOnly = Transcripts[-c(grep("ERCC-0", Transcripts$Feature)),]
   
   # Remove ERCCs in the definition file that are not in the count data file
   idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]
   
   # Remove ERCCs without a Ratio
   idCols = idCols[which(is.finite(idCols$Ratio)),]
   
   # Remove ERCCs from count data and idCols that are absent from the experiment
   TranscriptsERCCOnly = TranscriptsERCCOnly[match(idCols$Feature,
                                                  TranscriptsERCCOnly$Feature),]
   Transcripts = rbind(TranscriptsERCCOnly, TranscriptsHumanOnly)
   print(dim(Transcripts))
   print(length(grep("ERCC-",Transcripts$Feature)))
   ###############################################################################################################
   
   # Generate the Design Matrix for the table
   # example name is SEQC_ILM_BGI_A_1_L01_ATCACG_AC0AYTACXX
   designMat = getDesignMat(TranscriptsERCCOnly, 
                            factorList = c("Study","Platform","Site","Sample",
                                           "Library","Lane","Barcode",
                                           "Flowcell"), patternSplit = '_')
   if (sampleInfo$platform == "ROC"){
     designMat = getDesignMat(TranscriptsERCCOnly,
                              factorList = c("Sample","Library","Replicate",
                                             "Site"), patternSplit = '_')
   }
   # Subset just the A and B samples for Count Matrix and the designMatrix
   select = subset(designMat, (Sample == sampleInfo$sample1Name)|(Sample == sampleInfo$sample2Name))
   select <- as.data.frame(lapply(select,as.character))
   select <- as.data.frame(lapply(select,as.factor))
   designMatAll = select
   dataAB = Transcripts[-c(1)]
   TranscriptsAB = cbind(Transcripts[c(1)],dataAB[c(match(select$countSet,
                                                          names(dataAB)))])
   
   # Change the sample names to UHRR and HBRR first
   colnames(TranscriptsAB) <- gsub(pattern = "_A_",replacement="_UHRR_",
                                   x=colnames(TranscriptsAB))
   colnames(TranscriptsAB) <- gsub(pattern = "_B_",replacement="_HBRR_",
                                   x=colnames(TranscriptsAB))
   
   if(sampleInfo$platform == "ROC"){
     colnames(TranscriptsAB) <- gsub(pattern = "A_",replacement="UHRR_",
                                     x=colnames(TranscriptsAB))
     colnames(TranscriptsAB) <- gsub(pattern = "B_",replacement="HBRR_",
                                     x=colnames(TranscriptsAB))
     
   }
   sample1 = "UHRR";sample2 = "HBRR";
   
   designMatAB = getDesignMat(TranscriptsAB,
                              factorList = c("Study","Platform","Site","Sample",
                                             "Library","Lane","Barcode",
                                             "Flowcell"), patternSplit = '_')
   if (sampleInfo$platform == "ROC"){
     designMatAB = getDesignMat(TranscriptsAB,
                                factorList = c("Sample","Library","Replicate",
                                               "Site"), patternSplit = '_')
   }
   if(sampleInfo$platform != "ROC"){
     mainReads$sampleName <- gsub(pattern = "_A_", replacement="_UHRR_",
                                  x=mainReads$sampleName)
     mainReads$sampleName <- gsub(pattern = "_B_", replacement="_HBRR_",
                                  x=mainReads$sampleName)
     totalReads = mainReads$numReads[c(match(designMatAB$countSet,
                                             mainReads$sampleName))]  
   }else{
     mainReads$sampleName <- colnames(TranscriptsAB)[-1]
     totalReads = mainReads$TotalReads[c(match(designMatAB$countSet,
                                               mainReads$sampleName))]
     print(totalReads)
   }
   ### Reset TranscriptsAB and designMatAB to have the Lab# naming scheme
   colnames(TranscriptsAB) <- gsub(pattern = oldName, 
                                   replacement=sampleInfo$siteName,
                                   x=colnames(TranscriptsAB))
   designMatAB = getDesignMat(TranscriptsAB,
                              factorList = c("Study","Platform","Site","Sample",
                                             "Library","Lane","Barcode",
                                            "Flowcell"), patternSplit = '_')
   if (sampleInfo$platform == "ROC"){
     colnames(TranscriptsAB) <- gsub(pattern = oldName, 
                                     replacement=sampleInfo$siteName,
                                     x = colnames(TranscriptsAB))
     designMatAB = getDesignMat(TranscriptsAB,
                                factorList = c("Sample","Library","Replicate",
                                               "Site"), patternSplit = '_')
   }
   
   
 }
 if (sampleInfo$study == "SEQC_RatTox"){
   Transcripts <- read.delim("Data/RatToxSEQC/SEQC_TGx_GeneCounts_JMEEHAN.txt")
   # force names to be ERCC- and first column name to Feature
   names(Transcripts)[1] = "Feature"
   Transcripts$Feature = gsub("ERCC_","ERCC-",Transcripts$Feature)
   Transcripts$Feature = gsub(":Gene_AceView08","",Transcripts$Feature)
   Transcripts$Feature = gsub(":Gene_RefSeq","",Transcripts$Feature)
   # get the total reads per sample (from sequence files prior to mapping)
   ratToxReads <- read.csv("Data/RatToxSEQC/RatToxTotalReads54subset21.csv")
   
   #############################################################################
   
   # get data frames with just the ERCCs and just the human genes
   TranscriptsERCCOnly = Transcripts[c(grep("ERCC-0", Transcripts$Feature)),]
   TranscriptsHumanOnly = Transcripts[-c(grep("ERCC-0", Transcripts$Feature)),]
   
   # Remove ERCCs in the definition file that absent from count data file
   idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]
   
   # Remove ERCCs without a Ratio
   idCols = idCols[which(is.finite(idCols$Ratio)),]
   
   # Remove ERCCs from count data and idCols that are absent from the experiment
   TranscriptsERCCOnly = TranscriptsERCCOnly[match(idCols$Feature,
                                                  TranscriptsERCCOnly$Feature),]
   Transcripts = rbind(TranscriptsERCCOnly, TranscriptsHumanOnly)
   
################################################################################
   
   # Generate the Design Matrix for the table,
   # example name is SEQC_ILM_BGI_A_1_L01_ATCACG_AC0AYTACXX
   if(sampleInfo$analysis == "RatTox"){
     designMat = getDesignMat(TranscriptsERCCOnly, 
                              factorList = c("Site","PI","Flowcell","Barcode",
                                             "Tissue","Chemical","Vehicle",
                                             "Route","SeqBarcode","Lane"),
                              patternSplit = '_')
   }
   sample1 <- sampleInfo$sample1Name 
   sample2 <- sampleInfo$sample2Name
   # Subset to get A and B samples for the Count Matrix and for the designMatrix
   if ((sample1 == "NIT")|(sample1 == "THI")){
     select1 = subset(designMat,(Chemical == sample1)&(Vehicle == "NN")
                      &(Route == "IP"))
     print(select1)
     select2 = subset(designMat,(Chemical == sample2)&(Flowcell == "AB029JACXX")
                      &(Vehicle == "NN")&(Route == "IP")&
                        ((Lane == "s_4")|(Lane == "s_6")|(Lane == "s_1")))
     print(select2)
     select = rbind(select1,select2)  
   }
   if((sample1 == "3ME")|(sample1 == "MET")|(sample1 == "NAP")){
     select1 = subset(designMat,(Chemical == sample1)&(Vehicle == "NN")
                      &(Route == "OG"))
     print(select1)
     select2 = subset(designMat,(Chemical == sample2)&(Flowcell == "AB029JACXX")
                      &(Vehicle == "NN")&(Route == "OG"))
     print(select2)
     select = rbind(select1,select2)
   }
   
   select <- as.data.frame(lapply(select,as.character))
   select <- as.data.frame(lapply(select,as.factor))
   
   designMatAB = select
   
   dataAB = Transcripts[-c(1)]
   TranscriptsAB = cbind(Transcripts[c(1)],dataAB[c(match(select$countSet,
                                                          names(dataAB)))])
   totalReads = ratToxReads$total_reads[c(match(select$countSet,
                                                ratToxReads$Alt_ID))]
  
   dataAB <- TranscriptsAB[-c(1)]
   colnames(dataAB)<-paste(rep(c(sample1,sample2),each=ncol(dataAB)/2),
                           c(1:(ncol(dataAB)/2),1:(ncol(dataAB)/2)),sep="_")
   print(colnames(dataAB))
   
   TranscriptsAB <- cbind(Feature = TranscriptsAB$Feature, dataAB)
   designMatAll <- designMatAB
 }
   idxsample <- which((rowMeans(TranscriptsAB[-c(1)])>1)&(rowSums(
     TranscriptsAB[-c(1)]!=0)>=2))
   
   TranscriptsAB <- TranscriptsAB[idxsample,]
   
   TranscriptsAB$Feature <- as.factor(as.character(TranscriptsAB$Feature))
   
   measERCCs <- TranscriptsAB$Feature[grep("ERCC-0", TranscriptsAB$Feature)]
   
   insuffDat <- setdiff(idCols$Feature, measERCCs)
   
   print(paste("Transcripts were removed with a mean count < 1 or more than 2",
               "replicates with 0 counts.")) 
   print(paste("A total of",length(insuffDat),"out of",length(idCols$Feature),
               "ERCC controls were filtered"))
   print("The excluded ERCCs are:")
   print(insuffDat)
   print(paste("The remaining",length(measERCCs),"ERCC controls were analyzed"))   
  
  if(sampleInfo$analysis == "RatTox"){
    designMatSum <- getDesignMat(expressionData = TranscriptsAB,
                                 factorList = c("Sample","Replicate"),
                                 patternSplit = '_')
  }else{
    designMatSum <- NULL
  }
    
  # write Transcript csv file to directory
  #write.csv(TranscriptsAB, paste(sampleInfo$filenameRoot,"Transcripts.csv",sep="."),
  #          row.names = F)
  # collect everything to add to expDat
  expDat = append(expDat, list(TranscriptsAB = TranscriptsAB,
                               designMatAB = designMatAB, 
                               designMatSum = designMatSum,
                               sampleNames = c(sample1,sample2),
                               idCols = idCols,
                               totalReads = totalReads))
  return(expDat)
#   return(list(TranscriptsAB = TranscriptsAB, designMatAB = designMatAB, 
#                designMatAll = designMatAll, designMatSum = designMatSum,
#                sample1 = sample1, sample2 = sample2, idCols = idCols, 
#                totalReads = totalReads))

}