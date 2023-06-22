#!/bin/R
# By Anoushka Joglekar 04.2023
## Reads in multiway interaction data which contains cardinality of interactions,
## alignment coordinates, and closest gene to alignment. Also reads in a list of genes
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
readChromDF <- fread('../v0_poreC_explore/readChromunities_wClosestGene.gz')
colnames(readChromDF) <- c("ReadID","Cardinality","gene_ID","codingStatus","geneName","distToGene")

filteredReadChromDF <- readChromDF %>% filter(distToGene <= 100000)

## First get a null distribution -------
cardFreq <- readChromDF %>% select(ReadID,Cardinality,codingStatus) %>% 
  distinct() %>% 
  group_by(Cardinality,codingStatus) %>%
  summarise(numGenes = n())
  
ggplot(cardFreq %>% filter(codingStatus %in% c("lncRNA","protein_coding")), 
       aes(x = Cardinality, y = log10(numGenes), fill = codingStatus)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  xlim(0,30)

ggplot(readChromDF %>% filter(codingStatus %in% c("lncRNA","protein_coding")), 
       aes(x = Cardinality)) +
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlim(0,30)

## Protein coding only. Estimate the equation ----
install.packages('poweRlaw')
library(poweRlaw)
readChrom_pc <- readChromDF %>% filter(codingStatus == "protein_coding")
cardFreq_pc <- cardFreq %>% filter(codingStatus == "protein_coding")

ggplot(cardFreq_pc, 
       aes(x = Cardinality, y = numGenes)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlim(0,30)

plot(cardFreq_pc$Cardinality,cardFreq_pc$numGenes,xlim = c(0,30),pch = 20)
lines(spline(cardFreq_pc$Cardinality,cardFreq_pc$numGenes, n = 300, method = "natural"), col = 2)

f <- splinefun(cardFreq_pc$Cardinality,cardFreq_pc$numGenes, method = "natural")

counts <- cardFreq_pc$numGenes
pl <- displ$new(counts)
pl$setXmin(2)
pl$setPars(30)

plot(pl)
lines(pl)

## Trying to understand the data a bit more: --------
set.seed(1)
reads = sample(unique(readChrom_pc %>% filter(Cardinality == 10) %>% .$ReadID),100, replace = F)
df_subset <- readChromDF %>% filter(ReadID %in% reads)
se_subset <- left_join(df_subset,SE_df,by = "gene_ID")
ggplot(se_subset, aes(x = ReadID, y = Nuclear)) + 
  geom_point() + geom_line() + 
  theme(axis.text.x = element_blank())

## randomizing my cardinality assignments ---------
shuffledReadChromDF <- readChromDF
cards <- readChromDF %>% select(ReadID, Cardinality) %>% distinct()
cardIds <- cards$Cardinality
names(cardIds) <- cards$ReadID
shuffledReadChromDF$ReadID <- sample(shuffledReadChromDF$ReadID)
shuffledReadChromDF$Cardinality <- cardIds[shuffledReadChromDF$ReadID]

shufCombinedDF <- left_join(shuffledReadChromDF,SE_df,by = "gene_ID") %>% 
  group_by(ReadID) %>% mutate(medianSE=median(Nuclear,na.rm = T), 
                              missValues = sum(is.na(Nuclear))/n())

shufCombined_mod <- shufCombinedDF %>% dplyr::select(ReadID,Cardinality,codingStatus, medianSE,missValues) %>% distinct()

## Filter low cardinality and fewer missing values for splicing efficiency:
shufCombined_mod_filt <- shufCombined_mod %>% filter(!is.na(medianSE) & 
                                               Cardinality <= 10 & missValues <= 0.25)
shufCombined_mod_filt <- as.data.frame(shufCombined_mod_filt)
shufCombined_mod_filt$Cardinality <- as.factor(shufCombined_mod_filt$Cardinality)

dim(shufCombined_mod_filt)
ggplot(shufCombined_mod_filt, 
       aes(x = Cardinality, y = medianSE)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean splicing efficiency") +
  geom_signif(comparisons = list(c('2','3'),c('3','4'),c('4','5'),c('5','6'),c('6','7'),
                                 c('7','8'),c('8','9'),c('9','10')),
              map_signif_level = TRUE,tip_length = 1,
              y_position = c(1.04,1.03,1.02,1.01,1.04,1.03,1.02,1.01)) +
  ylim(0.75,1.1) 

## for all the possible cardinality ignoring missing values:

shufCombined_mod_filt <- shufCombined_mod %>% filter(!is.na(medianSE) & Cardinality < 30)
shufCombined_mod_filt <- as.data.frame(shufCombined_mod_filt)
shufCombined_mod_filt$Cardinality <- as.factor(shufCombined_mod_filt$Cardinality)

dim(shufCombined_mod_filt)
ggplot(shufCombined_mod_filt,
       aes(x = Cardinality, y = medianSE)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean splicing efficiency") +
  ylim(0.75,1) 

### Let's simulate a null expected distribution: ---------
## same number of points per cardinality:
null_sameDF <- do.call("rbind",lapply(1:13, function(i) do.call("rbind",lapply(1:100, function(s) 
  data.frame(sampleID = s, cardinality = i,
             meanSE = mean(sample(SE_df$Nuclear,i)))))))

null_sameDF$cardinality <- as.factor(null_sameDF$cardinality)
ggplot(null_sameDF,aes(x = cardinality, y = meanSE)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.75,1.1)

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

# Are there more interactions within prot-coding genes vs noncoding? ------
numberOfGenes <- readChromDF %>% select(geneName,codingStatus) %>% 
  filter(codingStatus %in% c("lncRNA","protein_coding")) %>% 
  distinct() %>% select(-geneName) %>% 
  group_by(codingStatus) %>% summarize(count = n())
ratio <- 0.637


bioType_interactions <- readChromDF %>% select(codingStatus,gene_ID,ReadID,Cardinality) %>%
  filter(codingStatus %in% c("protein_coding","lncRNA")) %>%
  select(-gene_ID) %>% group_by(ReadID, codingStatus) %>%
  add_count() %>% distinct() %>% ungroup() 

bioTypeIntSummary <- bioType_interactions %>% 
  group_by(codingStatus,Cardinality,n) %>%
       summarise(nCount = n()) %>% 
  ungroup() %>% group_by(Cardinality,codingStatus) %>% 
  mutate(Observed = nCount/sum(nCount)) %>%
  filter(Cardinality <= 20)

binomialProbs <- do.call('rbind',lapply(2:20, function(x) 
  data.frame(Expected = dbinom(1:x, x, ratio), Cardinality = x,n = 1:x)))

combinedBinomDF <- left_join(bioTypeIntSummary,binomialProbs) %>%
  pivot_longer(c("Expected","Observed"),names_to = "Type",values_to = "Percentage")

plotList <- NULL
for(i in 1:19){
  plotList[[i]] <- ggplot(combinedBinomDF %>% filter(Cardinality == (i+1) & codingStatus == "protein_coding"),
       aes(x = n, y = Percentage, col = Type)) +
  geom_point() + geom_line() + ylim(0,1) + 
    theme_minimal() + theme(legend.position = "none") +
    ggtitle(paste0("Cardinality: ",i+1)) +
    xlab("n prot coding genes") + ylab("% chromunities")
}

grid.arrange(grobs = plotList[1:15], ncol = 5)

## Correlate with splicing -------
## very basic QC 
combinedDF <- left_join(filteredReadChromDF,SE_df,by = "gene_ID") %>% 
  group_by(ReadID) %>% mutate(meanSE=mean(Nuclear,na.rm = T), 
                              missValues = sum(is.na(Nuclear))/n())

combined_mod <- combinedDF %>% dplyr::select(ReadID,Cardinality,codingStatus, meanSE,missValues) %>% distinct()

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



## for all the possible cardinality ignoring missing values:

combined_mod_filt <- combined_mod %>% filter(!is.na(meanSE) & Cardinality < 30)
combined_mod_filt <- as.data.frame(combined_mod_filt)
combined_mod_filt$Cardinality <- as.factor(combined_mod_filt$Cardinality)

dim(combined_mod_filt)
ggplot(combined_mod_filt,
       aes(x = Cardinality, y = meanSE)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean splicing efficiency") +
  ylim(0.75,1) 

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

## log odds ratio etc

logOddsDF_mean <- expAndObs %>% group_by(Cardinality,Type) %>%
  summarise(meanOfDistr = mean(meanSE)) %>%
  pivot_wider(names_from = Type, values_from = meanOfDistr) %>%
  mutate(OR_mean = Observed/Expected)
  
logOddsDF_median <- expAndObs %>% group_by(Cardinality,Type) %>%
  summarise(medianOfDistr = median(meanSE)) %>%
  pivot_wider(names_from = Type, values_from = medianOfDistr) %>%
  mutate(OR_median = Observed/Expected)

logOddsDF_mean$Cardinality <- as.integer(logOddsDF_mean$Cardinality)
ggplot(logOddsDF_mean, aes(Cardinality, OR_mean)) +
  geom_point() + geom_smooth() + ylim(1,1.2)

logOddsDF_median$Cardinality <- as.integer(logOddsDF_median$Cardinality)
ggplot(logOddsDF_median, aes(Cardinality, OR_median)) +
  geom_point() +
  geom_smooth() + ylim(1,1.2) +
  geom_text(aes(x = 10,y = 1.15), label = "p=0.003")

cor.test(logOddsDF_median$Expected,logOddsDF_median$Observed)
library(philentropy)

## Correlate with gene expression -------
gEx <- NULL
filePath <- '../2023_03_06_GM12878_cellularFractionData/GM12878_Nuclear_X/featureCountsOut/GM12878_Nuclear_X.gene.counts.txt'
files <- sapply(c("Rep1","Rep2"),function(rep) gsub("X", rep, filePath))
for (f in files){
  repName <- unlist(strsplit(basename(f),"_|[.]"))[3]
  gEx[[repName]] <- fread(f) %>% mutate(Rep = repName)
  colnames(gEx[[repName]])[1:2] <- c("gene_ID","counts")
}
gEx_df <- do.call('rbind',gEx) 

## plot gene expression concordance between replicates
plotGEX <- gEx_df %>% pivot_wider(names_from = Rep, values_from = counts)
ggplot(plotGEX,aes(x = log10(Rep1), y = log10(Rep2))) + 
  geom_point(size = 0.1, alpha = 0.5)

gEx_df <- gEx_df %>% group_by(Rep) %>% mutate(tot = sum(counts)) %>% 
  ungroup() %>%
  group_by(gene_ID) %>% summarize(meanCount = mean(counts*1e6/tot))

## Get mean expression of genes for read chromunities
combined_chromGEx_DF <- left_join(readChromDF,gEx_df,by = "gene_ID") %>% 
  group_by(ReadID) %>% mutate(meanExpressionPerRead=mean(log10(meanCount+0.5),na.rm = T))

## all
combined_chromGEx_mod <- combined_chromGEx_DF %>% select(ReadID,Cardinality, meanExpressionPerRead) %>% distinct()
combined_chromGEx_mod_filt <- combined_chromGEx_mod %>% filter(Cardinality < 30)
combined_chromGEx_mod_filt <- as.data.frame(combined_chromGEx_mod_filt)
combined_chromGEx_mod_filt$Cardinality <- as.factor(combined_chromGEx_mod_filt$Cardinality)

dim(combined_chromGEx_mod_filt)
ggplot(combined_chromGEx_mod_filt,
       aes(x = Cardinality, y = meanExpressionPerRead)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean (log) expression per chromunity") 

## Filter low cardinality:
combined_chromGEx_mod_filt <- combined_chromGEx_mod %>% filter(Cardinality <= 13)
combined_chromGEx_mod_filt <- as.data.frame(combined_chromGEx_mod_filt)
combined_chromGEx_mod_filt$Cardinality <- as.factor(combined_chromGEx_mod_filt$Cardinality)

dim(combined_chromGEx_mod_filt)
ggplot(combined_chromGEx_mod_filt,
       aes(x = Cardinality, y = meanExpressionPerRead)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("Mean (log) expression per chromunity") + 
  geom_signif(comparisons = list(c('2','3'),c('3','4'),c('4','5'),c('5','6'),c('6','7'),
                                 c('7','8'),c('8','9'),c('9','10')),
              map_signif_level = TRUE,tip_length = 0,
              y_position = c(5.05,5,4.05,4,5.05,5,4.05,4))


