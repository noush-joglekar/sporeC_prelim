#!/bin/R
# By Anoushka Joglekar 06.2023
## Reads in multiway interaction data which contains cardinality of interactions,
## alignment coordinates, closest gene to alignment, and read length. Also reads in a list of genes
## with splicing efficiency and expression counts per gene in the nuclear fraction.
## At present this is for GM12878 data only, but can be expanded to other cell lines / tissue

## This is basic QC to see if we see any positive / negative trends between cardinality and
## gene expression or splicing efficiency

library(data.table)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(tidyr)
library(gridExtra)

scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

SE_df <- read.table('../2023_03_06_GM12878_cellularFractionData/v0.evaluateSplicingEfficiency/SE_AvgOverReps_NuclearSE_withStatus', 
                    header = T) %>% dplyr::select(gene_ID,Nuclear,nucStatus)
readChromDF <- fread('../v0_poreC_explore/readChromunities_wClosestGeneAndReadLengths.gz')
colnames(readChromDF) <- c("ReadID","Cardinality","gene_ID","codingStatus","geneName","distToGene","readLength","readQual")

filteredReadChromDF <- readChromDF %>% filter(distToGene <= 100000)

### Let's simulate a null expected distribution: ---------
## different number of points per cardinality:
set.seed(10)
numPoints <- sort(sample(100:10000,12),decreasing = T)
null_diffDF <- do.call("rbind",lapply(1:12, function(i) do.call("rbind",lapply(1:numPoints[i], function(s) 
  data.frame(sampleID = s, cardinality = i+1,
             meanSE = mean(sample(SE_df$Nuclear,i+1)))))))

null_diffDF$cardinality <- as.factor(null_diffDF$cardinality)
ggplot(null_diffDF,aes(x = cardinality, y = meanSE)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.75,1.1)

bpStats_expected <- null_diffDF %>% group_by(cardinality) %>% 
  summarise(boxplot= list( setNames(boxplot.stats(meanSE)$stats,
                                    c('lower_whisker','lower_hinge','median','upper_hinge','upper_whisker')))) %>%
  unnest_wider(boxplot)


## Correlate with splicing -------
## very basic QC 
combinedDF <- left_join(filteredReadChromDF,SE_df,by = "gene_ID") %>% 
  group_by(ReadID) %>% mutate(meanSE=mean(Nuclear,na.rm = T), 
                              missValues = sum(is.na(Nuclear))/n())

combined_mod <- combinedDF %>% dplyr::select(ReadID,Cardinality,codingStatus,readLength,readQual,
                                             meanSE,missValues) %>% distinct()

## Filter low cardinality and fewer missing values for splicing efficiency:
combined_mod_filt <- combined_mod %>% filter(!is.na(meanSE) & 
                                               Cardinality <= 13 & missValues <= 0.25)
combined_mod_filt <- as.data.frame(combined_mod_filt)
combined_mod_filt$Cardinality <- as.factor(combined_mod_filt$Cardinality)

dim(combined_mod_filt)
ggplot(combined_mod_filt, 
       aes(x = Cardinality, y = meanSE)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean splicing efficiency") +
  geom_signif(comparisons = list(c('2','3'),c('3','4'),c('4','5'),c('5','6'),c('6','7'),
                                 c('7','8'),c('8','9'),c('9','10')),
              map_signif_level = TRUE,tip_length = 1,
              y_position = c(1.04,1.03,1.02,1.01,1.04,1.03,1.02,1.01)) +
  ylim(0.75,1.1)

ggplot(combined_mod_filt, 
       aes(x = Cardinality, y = meanSE)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean splicing efficiency") +
  ylim(0.75,1)


## Boxplot of cardinality vs. read length
ggplot(combined_mod_filt, 
       aes(x = Cardinality, y = readLength)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Read length") +
  ylim(0,15000)

## get boxplot stats
table(combined_mod_filt$Cardinality)
bpStats <- combined_mod_filt %>% group_by(Cardinality) %>% 
  summarise(boxplot= list( setNames(boxplot.stats(meanSE)$stats,
                                    c('lower_whisker','lower_hinge','median','upper_hinge','upper_whisker')))) %>%
  unnest_wider(boxplot)

obs_overE <- data.frame(Cardinality = 2:13,ratioOfMedians = bpStats$median/bpStats_expected$median)


## make a decent-ish figure panel
expAndObs <- rbind(combined_mod_filt %>% dplyr::select(ReadID, Cardinality, meanSE) %>% mutate(Type = "Observed"), 
                   null_diffDF %>% rename(ReadID = sampleID, Cardinality = cardinality) %>% mutate(Type = "Expected"))
g0 <- ggplot(SE_df, aes(x = Nuclear)) + 
  geom_histogram(fill = "#FEE6CE") +
  theme_minimal(base_size = 16) + xlab("Splicing efficiency") 

g1 = ggplot(expAndObs, 
            aes(x = Cardinality, y = meanSE, fill = Type)) + 
  geom_boxplot(notch = T, outlier.shape = NA, coef = 0.5) + 
  ylab("Mean splicing efficiency") +
  scale_fill_brewer(palette = 7) + 
  ylim(0.75,1) + theme_minimal(base_size = 16) +
  theme(legend.position = "bottom")

g2 <- ggplot(obs_overE, aes(x = Cardinality, y = ratioOfMedians)) +
  geom_point(col = "#FDAE6B")+ theme_minimal(base_size = 14) +
  ylim(1,1.2) + geom_smooth(col = "#FDAE6B") + 
  geom_text(aes(x = 10,y = 1.18), label = "p=0.003\nPearson rho =0.77") +
  ylab("odds Ratio of Medians")

g1 / (g0 | g2)