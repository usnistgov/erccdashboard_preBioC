### Evaluate the sd(ratios) for the ERCCs and the simulated genes
# Hypothesis is that there is a difference between the empirical and simulated
# results

## this needs work ...
simcnt <- expDat$simcnt
sdSimPlot <- NULL
for (i in 1:nrow(simcnt)){
  sdSim <- sd(log2(simcnt[i, 1:4]/simcnt[i,5:8]))
  meanSim <- mean(log2(simcnt[i,]))
  row <- data.frame(sdSim, meanSim)
  sdSimPlot <- rbind(sdSimPlot, row)
}
ggplot(sdSimPlot) + geom_point(aes(x = meanSim, y = sdSim))