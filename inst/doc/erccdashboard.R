### R code from vignette source 'erccdashboard.Rnw'

###################################################
### code chunk number 1: erccdashboard.Rnw:14-18
###################################################
options(width=60, continue = "  ")
#options(SweaveHooks=list(fig=function()
#                par(mar=c(5.1,4.1,1.1,2.1))))
library( "erccdashboard" )


###################################################
### code chunk number 2: loadExampleData
###################################################
data(SEQC.Example)


###################################################
### code chunk number 3: defineInputData
###################################################
datType = "count" # "count" for RNA-Seq data, "array" for microarray data
isNorm = FALSE # flag to indicate if input expression measures are already
               # normalized, default is FALSE 
exTable = COH.RatTox.ILM.MET.CTL.countTable # the expression measure table
filenameRoot = "COH.ILM" # user defined filename prefix for results files
sample1Name = "MET" # name for sample 1 in the experiment
sample2Name = "CTL" # name for sample 2 in the experiment
erccmix = "RatioPair" # name of ERCC mixture design, "RatioPair" is default
erccdilution = 1/100 # dilution factor used for Ambion spike-in mixtures
spikeVol = 1 # volume (in microliters) of diluted spike-in mixture added to 
             #   total RNA mass
totalRNAmass = 0.500 # mass (in micrograms) of total RNA 
choseFDR = 0.05 # user defined false discovery rate (FDR), default is 0.05


###################################################
### code chunk number 4: inspectRatCount
###################################################
head(COH.RatTox.ILM.MET.CTL.countTable)


###################################################
### code chunk number 5: runDashboardSEQCcount
###################################################
exDat <- runDashboard(datType = "count", isNorm = FALSE,
                       exTable = COH.RatTox.ILM.MET.CTL.countTable,
                       filenameRoot = "COH.ILM",sample1Name = "MET",
                       sample2Name = "CTL",erccmix = "RatioPair",
                       erccdilution = 1/100, spikeVol = 1,
                       totalRNAmass = 0.500, choseFDR = 0.1)


###################################################
### code chunk number 6: initializeData
###################################################
summary(exDat)


###################################################
### code chunk number 7: ratPlotA
###################################################
exDat$Figures$dynRangePlot


###################################################
### code chunk number 8: ratPlotB
###################################################
exDat$Figures$rocPlot


###################################################
### code chunk number 9: ratPlotC
###################################################
exDat$Figures$lodrERCCPlot


###################################################
### code chunk number 10: ratPlotD
###################################################
exDat$Figures$maPlot


###################################################
### code chunk number 11: erccdashboard.Rnw:237-243
###################################################
exDat <- runDashboard(datType="count", isNorm = FALSE,
                      exTable=Lab5.ILM.UHRR.HBRR.countTable,
                      filenameRoot="Lab5", sample1Name = "UHRR", 
                      sample2Name = "HBRR", erccmix = "RatioPair",
                      erccdilution = 1, spikeVol = 50, 
                      totalRNAmass = 2.5*10^(3), choseFDR = 0.01)


###################################################
### code chunk number 12: seqcCountPlotA
###################################################
exDat$Figures$dynRangePlot


###################################################
### code chunk number 13: totalReadsCompare
###################################################
COH.RatTox.ILM.MET.CTL.totalReads
Lab5.ILM.UHRR.HBRR.totalReads


###################################################
### code chunk number 14: seqcCountPlotB
###################################################
exDat$Figures$rocPlot


###################################################
### code chunk number 15: seqcCountPlotC
###################################################
exDat$Figures$lodrERCCPlot


###################################################
### code chunk number 16: seqcCountPlotD
###################################################
exDat$Figures$maPlot


###################################################
### code chunk number 17: SEQCMainArray
###################################################
exDat <- runDashboard(datType="array", isNorm = F,
                      exTable=Lab13.array.UHRR.HBRR,
                      filenameRoot = "Lab13.array",
                      sample1Name = "UHRR", sample2Name="HBRR",
                      erccmix = "RatioPair", erccdilution = 1, 
                      spikeVol = 50, totalRNAmass = 2.5*10^(3), choseFDR=0.01)


###################################################
### code chunk number 18: seqcArrayPlotA
###################################################
exDat$Figures$dynRangePlot


###################################################
### code chunk number 19: seqcArrayPlotB
###################################################
exDat$Figures$rocPlot


###################################################
### code chunk number 20: seqcArrayPlotC
###################################################
exDat$Figures$lodrERCCPlot


###################################################
### code chunk number 21: seqcArrayPlotD
###################################################
exDat$Figures$maPlot


###################################################
### code chunk number 22: viewDashboardOrder
###################################################
runDashboard


###################################################
### code chunk number 23: sessionData
###################################################
sessionInfo()


