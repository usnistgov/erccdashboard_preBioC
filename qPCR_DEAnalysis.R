### Analyze Taqman data

# load the ABRF processed version of the MAQC data. Checked with orignal MAQC GPL4097 family.soft files to make sure that this table from ABRF has been normalized to POLR2A

ABRFtaqman = read.table(file="../ERCCDashboardDevelopment/Manuscript/131219NBT_ThirdTime/ABRFdata/TaqmanPCR.txt",header=T)

colnames(ABRFtaqman)<-gsub("[.]","_",colnames(ABRFtaqman))

y <- ABRFtaqman
array1 <- as.character(y$gene_symbol)
array2 <- as.character(seq(from=1,to=length(array1),by=1))

row.names(y) <- paste(array1,array2,sep="_")

y <- y[-c(1)]

y <- y[c(1:8)]
design <- cbind(Grp1=1,Grp1vs2=c(1,1,1,1,0,0,0,0))

fit <- lmFit(y,design)
###fit <- lmFit(ylog,design)

fit <- eBayes(fit)

res <- topTable(fit,sort.by="none",number = dim(y)[1],coef = 2)

# use qvalue package to create qvals from P.vals
library(qvalue)
pval <- res$adj.P.Val
qobj <- qvalue(pval)
res$qvals <- qobj$qvalues

# This is the FDR used for all MAQC sites
choseFDR <- 0.01

if(any(res$qvals<choseFDR)){
  p.thresh<-max(res$adj.P.Val[res$qvals<choseFDR])
}

# Alternative option is to set p.thresh directly based on work in MAQC
# e.g. choice of p = ??
res.DE <- subset(res, (abs(logFC)>=1) & (adj.P.Val <= p.thresh) )