#!/bin/R 
# By Anoushka Joglekar 08.2023
## Take synergistic interactions and evaluate splicing efficiency in those hubs 

## Setup -----
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/'
setwd(workingDir)

## Load libraries -------
library(chromunity)
library(gUtils)
library(data.table)
library(dplyr)
library(ggplot2)
library(tibble)
library(edgeR)

## Read in input -----
args <- commandArgs(trailing = TRUE)
args <- c('synergyRun_byChrom/allSynergies_chr12.tab','synergyRun_byChrom/slidingWinChrom_chr12.rds')
synEvents <- fread(args[1])
allChromunities <- readRDS(args[2])
SE_genes <- read.table('../v1_poreC_explore/Gene_NucSE',sep = "\t")
colnames(SE_genes) <- c("GeneID","Gene","SE")

rep1 <- read.table('../2023_03_06_GM12878_cellularFractionData/GM12878_Nuclear_Rep1/featureCountsOut/GM12878_Nuclear_Rep1.gene.counts.txt',
                   header = TRUE)
rep2 <- read.table('../2023_03_06_GM12878_cellularFractionData/GM12878_Nuclear_Rep2/featureCountsOut/GM12878_Nuclear_Rep2.gene.counts.txt',
                   header = TRUE)
gEx <- left_join(rep1,rep2) %>% column_to_rownames("Geneid")
colnames(gEx) <- c("Rep1","Rep2")
norm_gEx <- cpm(gEx)

norm_gEx_df <- as.data.frame(norm_gEx) %>% rownames_to_column("GeneID")


## Start processing ------
set.seed(10)
sigSyn <- synEvents %>% filter(fdr <= 0.1)
nonsigSyn <- synEvents %>% filter(fdr > 0.8) %>% slice_sample(n = nrow(sigSyn))

sigConcatemers <- gr2dt(allChromunities$concatemers %Q% (chid %in% sigSyn$bid)) %>%
  mutate(Status = "Sig")

# nonSigConcatemers <- gr2dt(allChromunities$concatemers %Q% (chid %in% nonsigSyn$bid)) %>%
#   mutate(Status = "NonSig")
untestedEvents = sample(setdiff(unique(allChromunities$concatemers$chid),synEvents$bid),size = nrow(sigSyn))
nonSigConcatemers <- gr2dt(allChromunities$concatemers %Q% (chid %in% untestedEvents)) %>%
  mutate(Status = "NonSig")

concatemerSet <- rbind(sigConcatemers,nonSigConcatemers)

## run some stats
numSig <- nrow(sigSyn)
numReadsPerSyn <- sigConcatemers %>% select(read_idx,chid) %>% distinct() %>% group_by(chid) %>% summarize(numReads = n())
avgCardPerSyn <- sigConcatemers %>% select(count,chid) %>% distinct() %>% group_by(chid) %>% summarize(avgCard = mean(count))
avgGenesPerSyn <- sigConcatemers %>% select(Gene,read_idx,chid) %>% distinct() %>% group_by(chid,read_idx) %>% add_count() %>% 
  ungroup() %>% select(-read_idx) %>% group_by(chid) %>% summarize(avgGenes = mean(n))

plot(avgCardPerSyn$avgCard,avgGenesPerSyn$avgGenes, pch = 16, xlim = c(2,18), ylim = c(2,18))

## SE per synergy
concatemerSet_wSE <- left_join(concatemerSet,SE_genes,by = "Gene") %>% tidyr::drop_na()
concatemerSet_wSE_gEx <- left_join(concatemerSet_wSE,norm_gEx_df,by = "GeneID") %>% 
  rowwise() %>% mutate(logCPM = log10(mean(c(Rep1,Rep2))))

meanSE_perRead <- concatemerSet_wSE_gEx %>% group_by(chid,read_idx) %>% 
  mutate(avgSE = mean(SE,na.rm = T), 
         dispSE = max(max(SE) - median(SE), median(SE)-min(SE)),
         avgGEx = mean(logCPM,na.rm = T)) %>% 
  select(chid,read_idx,count,avgSE,dispSE,avgGEx,Status) %>% distinct()

ggplot(meanSE_perRead %>% filter(count <= 15), 
       aes(x = as.factor(count),y = avgSE, fill = Status)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.75,1) 

ggplot(meanSE_perRead %>% filter(count <= 15), 
       aes(x = as.factor(count),y = dispSE, fill = Status)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal(base_size = 18)

ggplot(meanSE_perRead %>% filter(count <= 15), 
       aes(x = as.factor(count),y = avgGEx, fill = Status)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal(base_size = 18)

numReadsPerCard = concatemerSet_wSE_gEx %>% filter(Status == "Sig") %>% 
  select(count,chid) %>% group_by(count) %>% add_count() %>% 
  select(-chid) %>% distinct() %>% filter(count <= 15)

## getSE_df_perCard
combinedDF_chr19 <- combinedDF %>% select(ReadID,Cardinality,meanSE) %>% mutate(Status = "Background")
meanSE_redefined <- meanSE_perRead %>% rename(ReadID = read_idx, Cardinality = count, meanSE = avgSE) %>% 
  filter(Status == "Sig") %>% ungroup() %>% select(ReadID,Cardinality,meanSE,Status)
sigAndBG <- rbind(meanSE_redefined,combinedDF_chr19)

ggplot(sigAndBG %>% filter(Cardinality >=3 & Cardinality <= 15), 
       aes(x = as.factor(Cardinality),y = meanSE, fill = Status)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal(base_size = 18) +
  ylim(0.75,1)

ggplot(sigAndBG %>% filter(Cardinality <= 15), 
       aes(x = Status,y = meanSE, fill = Status)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal(base_size = 18) +
  ylim(0.75,1)






