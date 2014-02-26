#microarray test commands

load("data/SEQC.Prelim.Array.RData")
expTable <- UTSW.array.UHRR.HBRR
repNormFactor <- NULL
filenameRoot = "UTSW.Array"
sample1Name = "UHRR"
sample2Name = "HBRR"
ERCCmixes = "RatioPair"
ERCCdilution = 1/100
spikeVol = 1
totalRNAmass = 0.500
choseFDR = 0.01
datType = "array"

expDat <- initDat(datType, expTable, repNormFactor, filenameRoot, sample1Name,
                  sample2Name, ERCCmixes, ERCCdilution, spikeVol, totalRNAmass,
                  choseFDR)

expDat <- est_r_m(expDat)

#set up different p-value function
expDat <- geneExprTest(expDat)

#
expDat <- erccROC(expDat)

#
expDat = estLODR(expDat,kind = "ERCC", prob=0.9)

#
expDat <- dynRangePlot(expDat, errorBars = T)

#
expDat <- annotLODR(expDat)

save(expDat,file=paste0(expDat$sampleInfo$filenameRoot,".RData"))

savePlots(expDat)

