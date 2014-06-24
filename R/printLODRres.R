printLODRres <- function(exDat){
    sampleInfo <- exDat$sampleInfo
    erccInfo <- exDat$erccInfo
    
    fit.coeff <- exDat$fit.coeff
    mnLibeFactor <- exDat$mnLibeFactor
    
    FCcode = erccInfo$FCcode
    legendLabels = erccInfo$legendLabels
    
    
    lodr.res = data.frame(exDat$Results$lodr.res.ERCC)
    #print(lodr.res)
    if(dim(lodr.res)[2]!=0){
        ### Fold = lodr.res[c(1)]
        Fold = as.numeric(exDat$erccInfo$FCcode$FC)
        
        Abund = as.numeric(gsub("<", "",lodr.res$Estimate))
        Ratio = legendLabels
        
        if(is.null(mnLibeFactor)){
            # Convert LODR Abund estimate to library size normalized
            logAbund = log2(Abund)
        }else{
            # Convert LODR Abund estimate to library size normalized
            logAbund = log2((Abund/(mnLibeFactor)))
        }
        
        ###LODR.print.res = data.frame(Fold, Ratio, Abund, logAbund, ConcEst)
        LODR.print.res = data.frame(Fold, Ratio, Abund, logAbund)
        
        names(LODR.print.res)<- c("Fold","Ratio","Abund","Log2Abund")
        #print(LODR.print.res)
        
        #cutoffs = ConcEst[which(!(is.na(ConcEst)))]
        AbundCutoffs = logAbund[which(!(is.na(logAbund)))]
    }else{
        LODR.print.res <- NULL
        AbundCutoffs <- NULL
    }
    
    return(list(LODRtable = LODR.print.res,
                #cutoffs=cutoffs,
                AbundCutoffs = AbundCutoffs))
    
}