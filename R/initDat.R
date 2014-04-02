#' Initialize the expDat list
#'
#' @param datType       type is "count" or "array", unnormalized data is  
#'                      expected (normalized data may be accepted in future
#'                      version of the package). Default is "count" (integer 
#'                      count data),"array" is unnormalized fluorescent 
#'                      intensities from microarray
#'                      fluorescent intensities (not log transformed or 
#'                      normalized)
#' @param expTable      data frame, the first column contains names of 
#'                      genes or transcripts (Feature) and the remaining columns
#'                      are counts for sample replicates spiked with ERCC 
#'                      controls
#' @param repNormFactor vector of normalization factors for each replicate
#' @param filenameRoot  string root name for output files
#' @param sample1Name   string name for sample 1 in the gene expression 
#'                      experiment
#' @param sample2Name   string name for sample 2 in the gene expression
#'                      experiment
#' @param erccmix     Name of ERCC mixture design, "RatioPair" is 
#'                      default, the other option is "Single"
#' @param erccdilution  unitless dilution factor used in dilution of the Ambion 
#'                      ERCC spike-in mixture solutions 
#' @param spikeVol      volume in microliters of diluted ERCC mix spiked into
#'                      the total RNA samples
#' @param totalRNAmass  mass in micrograms of total RNA spiked with diluted ERCC
#'                      mixtures 
#' @param choseFDR      False Discovery Rate for differential expression testing
#'                      , default is 0.05
#' @param ratioLim      Limits for ratio axis on MA plot, default is c(-4,4)
#' @param signalLim      Limits for ratio axis on MA plot, default is c(-12,12)
#' @param userMixFile   optional filename input, default is NULL, if ERCC 
#'                      control ratio mixtures other than the Ambion product
#'                      were used then a userMixFile can be used for the analysis
#'                                         
#' 
#' @export


initDat <- function(datType=NULL, expTable=NULL, repNormFactor=NULL,
                    filenameRoot = NULL,
                    sample1Name = NULL,sample2Name = NULL, 
                    erccmix = "RatioPair", erccdilution = 1,
                    spikeVol = 1, totalRNAmass = 1,choseFDR = 0.05,
                    ratioLim = c(-4,4), signalLim = c(-14,14), 
                    userMixFile =NULL){
  cat("\nInitializing the expDat list structure...\n")
  
  myYLimMA <- ratioLim
  myXLimMA <- signalLim
  
  myYLim <- myXLimMA
  myXLim <- NULL
  
  expDat<-NULL
  
  #myXLimMA = c(-10,15)
  #myYLimMA = c(-4,4)
#   
#   if ((datType == "count")|(datType == "array")){
#     myXLim = c(-10,15)
#   }else{
#     stop("datType is not count or array")
#     #need to define what x-axis will look like for FPKM
#     #chooseXLim <- function(){
#     #  cat("\nChoose X-scale, e.g. c(2,10)\n")
#     #  readline("Enter X-scale vector: ")
#     #}
#     #myXLim = as.numeric(chooseXLim())  
#   }
#   
  #myYLim = myXLimMA
  
  cat(paste("choseFDR =",choseFDR,"\n"))

  if(missing(userMixFile)){
    userMixFile <- NULL
  }
#   if((datType == "count") & (is.null(repNormFactor))){
#     stop("repNormFactor argument is missing!")
#   }
  if(is.null(repNormFactor)){
    #repNormFactor <- NULL
    cat("repNormFactor is NULL \n")
  }
  ## Do some library loading
  #library("QuasiSeq")
  #library("ROCR")
  #library("edgeR")
  #library("grid")
  #library("gridExtra")
  #library("reshape2")
  #library("gtools")
  #library("qvalue")
  
  ##############################
  libeSizeNorm = T ## set default to library size normalization
  sampleInfo = list(sample1Name = sample1Name,
                    sample2Name = sample2Name, choseFDR = choseFDR,
                    erccdilution = erccdilution, erccmix = erccmix,
                    spikeVol = spikeVol, totalRNAmass = totalRNAmass,
                    libeSizeNorm = libeSizeNorm, datType = datType)
  
  plotInfo = list(myXLimMA = myXLimMA, myYLimMA = myYLimMA,
                    myXLim = myXLim, myYLim = myYLim, xlimEffects = NULL)

 
  expDat <- list(sampleInfo = sampleInfo,plotInfo = plotInfo)
  
  if (exists("filenameRoot")){
    expDat <- dashboardFile(expDat,filenameRoot = filenameRoot)  
  }else{
    stop("The filenameRoot character string has not been defined!")
  }
  
  
  ###############################################################################
  # Run loadERCCInfo function to obtain ERCC information
  expDat <- loadERCCInfo(expDat, erccmix, userMixFile)
  
  ###############################################################################
  # Assume user has created data frame countTable and totalReads vector
  # process those data files to add to expDat structure
  expDat <- loadExpMeas(expDat, expTable, repNormFactor)
  
  

  ###############################################################################
  # library size normalize the data
  
  if(libeSizeNorm == T){
    expDat <- libeSizeNorm(expDat)  
  }
  
  
  ###############################################################################
  # length normalize the ERCC concentrations
  expDat <- prepERCCDat(expDat)
  
  # Estimate the mean library size factor for the data to use to estimate
  # corresponding concentrations for LODR
  #expDat <- estMnLibeFactor(expDat, cnt = expDat$Transcripts)
  
  expDat <- plotAdjust(expDat)
  
  return(expDat)
  
}