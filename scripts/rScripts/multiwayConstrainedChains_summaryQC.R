#!/bin/R

## Doing a quick QC on the chains with imposed multiway contacts
## because something is weird

workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/v3.multiwayConstraints/'
setwd(workingDir)

c1_summ <- as.data.frame(data.table::fread('cell2_output/summary.txt',keepLeadingZeros = T))
colnames(c1_summ) <- c("chainID","numReads","maxCard","primCutoff","secCutoff","diagDist","multiCode","cellType")

library(ggplot2)

hist(c1_summ$maxCard)
hist(log10(c1_summ$numReads))

mS = as.data.frame(table(c1_summ$multiCode))
colnames(mS)[1] <- "multiwayStatus"
ggplot(mS, aes(x = multiwayStatus, y = Freq)) +
  geom_bar(stat = "identity")

## Definitely weird. Going back to look at the distance matrices

## Now getting an idea of what the distance distribution looks like
library(tidyr)
library(dplyr)

cell1_summary <- read.csv('cell1_summaryDF_maxDists.csv') 

cell1_subDF <- as.data.frame(cell1_summary) %>% 
  dplyr::select(contains("Cell1")) %>% 
  pivot_longer(cols = contains("Inter"),names_to = "Interaction",values_to = "Dist") %>%
  mutate(Interaction = gsub("Cell1_","",Interaction), CellID = "Cell1") %>%
  rename(Status = Cell1_Status)

cell2_subDF <- as.data.frame(cell1_summary) %>% 
  dplyr::select(contains("Cell2")) %>% 
  pivot_longer(cols = contains("Inter"),names_to = "Interaction",values_to = "Dist") %>%
  mutate(Interaction = gsub("Cell2_","",Interaction), CellID = "Cell2") %>%
  rename(Status = Cell2_Status)

cell1_final <- rbind(cell1_subDF,cell2_subDF)
cell1_final$Status <- as.factor(cell1_final$Status)

ggplot(cell1_final, aes(x = Dist, fill = interaction(Interaction,CellID) )) +
  geom_histogram(alpha = 0.6,bins = 20) + 
  facet_wrap(~CellID) +
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom")

ggplot(cell1_final, aes(x = Dist, fill = CellID )) +
  geom_histogram(alpha = 0.6,bins = 20) + 
  theme_classic(base_size = 15) +
  theme(legend.position = "bottom")

## Look at joint distributions of distances 

cutoff = 300

c1_summDF = cell1_summary %>% 
  rowwise() %>% 
  mutate(c1m = min(c_across(contains("Cell1_Inter"))),
         c2m = min(c_across(contains("Cell2_Inter")))) %>%
  select(c1m,c2m) %>%
  mutate(status = case_when(c1m <= cutoff & c2m <= cutoff ~ '11',
                            c1m <= cutoff & c2m > cutoff ~ '10',
                            c1m > cutoff & c2m <= cutoff ~ '01',
                            c1m > cutoff & c2m > cutoff ~ '00'))


ggplot(c1_summDF, aes(x = c1m - c2m, fill = status)) +
  geom_histogram(alpha = 0.6) 

ggplot(c1_summDF, aes(x = c1m, fill = status)) +
  geom_histogram(alpha = 0.6) 
