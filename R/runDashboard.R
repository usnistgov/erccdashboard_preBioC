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
  # Required for all subsequent functions
  expDat <- initDat(datType=datType, expTable=expTable, 
                    repNormFactor=repNormFactor, filenameRoot=filenameRoot,
                    sample1Name=sample1Name, sample2Name=sample2Name, 
                    erccmix=erccmix, erccdilution=erccdilution, 
                    spikeVol=spikeVol, totalRNAmass=totalRNAmass,
                    choseFDR=choseFDR,userMixFile=userMixFile)
  
  # Estimate the difference in mRNA fraction of total RNA for the two samples
  # Required for all subsequent functions
  expDat <- est_r_m(expDat)
  
  # Evaluate the dynamic range of the experiment (Signal-Abundance plot)
  # Not required for subsequent functions
  expDat <- dynRangePlot(expDat)
  
  # Test for differential expression between samples
  # Required for all subsequent functions
  expDat <- geneExprTest(expDat)
  
  # Generate ROC curves and AUC statistics
  # Not Required for subsequent functions
  expDat <- erccROC(expDat)
  
  # Estimate LODR for ERCC controls
  # Required for subsequent functions
  expDat = estLODR(expDat,kind = "ERCC", prob=0.9)
  
  # Estimate LODR using Simulated data from endogenous transcripts
  # Not required for subsequent functions
  expDat = estLODR(expDat,kind = "Sim", prob=0.9)
  
  # Generate MA plot (Ratio vs. Average Signal) with ERCC controls below LODR 
  #   annotated also flags possible False Negatives on DE gene list based on LODR 
  #   threshold from DE gene list
  # Not required for subsequent functions
  expDat <- annotLODR(expDat)
  
  
  ### Saving plots and results
  # Convenience function to save 4 main figures to PDF
  saveERCCPlots(expDat)
  
  # Save expDat to a RData file for later use
  cat("\nSaving expDat list to .RData file...")
  nam <- paste(expDat$sampleInfo$filenameRoot, "expDat",sep = ".")
  assign(nam,expDat)
  
  to.save <- ls()
  
  save(list = to.save[grepl(pattern = nam,x=to.save)],
       file=paste0(expDat$sampleInfo$filenameRoot,".RData"))
  
  # End analysis and return expDat to global environ. / workspace
  cat("\nAnalysis completed.")
  return(expDat)
  
}