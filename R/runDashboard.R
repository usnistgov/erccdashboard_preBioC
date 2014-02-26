#' Run default erccdashboard analysis of ERCC control ratio mixtures
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
#' @param userMixFile   optional filename input, default is NULL, if ERCC 
#'                      control ratio mixtures other than the Ambion product
#'                      were used then a userMixFile can be used for the analysis
#'                                         
#' 
#' @export
#' @examples
#' load(file = system.file("data/SEQC.Example.RData",
#'      package = "erccdashboard"))
#'      
#' expDat = runDashboard(datType = "count",
#'                  expTable = COH.RatTox.ILM.MET.CTL.countTable, 
#'                  repNormFactor = COH.RatTox.ILM.MET.CTL.totalReads, 
#'                  filenameRoot = "COH.ILM",
#'                  sample1Name = "MET", sample2Name = "CTL", 
#'                  erccmix = "RatioPair", erccdilution = 1/100, 
#'                  spikeVol = 1, totalRNAmass = 0.500,choseFDR = 0.1)
#'                  
#' summary(expDat)

runDashboard <- function(datType=NULL, expTable=NULL, repNormFactor=NULL,
                         filenameRoot = NULL,
                         sample1Name = NULL,sample2Name = NULL, 
                         erccmix = "RatioPair", erccdilution = 1,
                         spikeVol = 1, totalRNAmass = 1,choseFDR = 0.05,
                         userMixFile=NULL){

  # Initialize expDat structure
  expDat <- initDat(datType=datType, expTable=expTable, 
                    repNormFactor=repNormFactor, filenameRoot=filenameRoot,
                    sample1Name=sample1Name, sample2Name=sample2Name, 
                    erccmix=erccmix, erccdilution=erccdilution, 
                    spikeVol=spikeVol, totalRNAmass=totalRNAmass,
                    choseFDR=choseFDR,userMixFile=userMixFile)
  
  # Estimate mRNA fraction of total RNA difference for two samples
  expDat <- est_r_m(expDat)
  
  # Test for differential expression
  expDat <- geneExprTest(expDat)
  
  # Generate ROC curves and AUC statistics
  expDat <- erccROC(expDat)
  
  # Estimate LODR for ERCC controls
  expDat = estLODR(expDat,kind = "ERCC", prob=0.9)
  
  # Estimate LODR using Simulated data from endogenous transcripts
  expDat = estLODR(expDat,kind = "Sim", prob=0.9)
  
  # Evaluate dynamic range of experiment with Signal-Abundance plot
  expDat <- dynRangePlot(expDat, errorBars = T)
  
  # Generate MA plot (Ratio vs. Average Signal) with ERCC controls below LODR 
  #   annotated also flag genes based on LODR threshold from DE gene list
  expDat <- annotLODR(expDat)
  
  
  #
  saveERCCPlots(expDat)
  
  cat("\nSaving expDat list to .RData file...")
  nam <- paste(expDat$sampleInfo$filenameRoot, "expDat",sep = ".")
  assign(nam,expDat)
  
  to.save <- ls()
  
  save(list = to.save[grepl(pattern = nam,x=to.save)],
       file=paste0(expDat$sampleInfo$filenameRoot,".RData"))
  
  
  cat("\nAnalysis completed.")
  return(expDat)
  
}