#!/bin/R

## Doing a quick QC on the chains with imposed multiway contacts
## because something is weird

workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/v3.multiwayConstraints/'
setwd(workingDir)

c1_summ <- as.data.frame(data.table::fread('cell1_output/summary.txt',keepLeadingZeros = T))
colnames(c1_summ) <- c("chainID","numReads","maxCard","primCutoff","secCutoff","diagDist","multiCode","cellType")

library(ggplot2)

hist(c1_summ$maxCard)
hist(log10(c1_summ$numReads))

mS = as.data.frame(table(c1_summ$multiCode))
colnames(mS)[1] <- "multiwayStatus"
ggplot(mS, aes(x = multiwayStatus, y = Freq)) +
  geom_bar(stat = "identity")

## Definitely weird. Going back to look at the distance matrices