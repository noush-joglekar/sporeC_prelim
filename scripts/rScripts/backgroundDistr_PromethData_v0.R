# Library loading ------
library(data.table)
library(dplyr)
library(tidyr)

## Setup
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/'
setwd(workingDir)


## Read in input
SE_df <- read.table('../2023_03_06_GM12878_cellularFractionData/v0.evaluateSplicingEfficiency/SE_AvgOverReps_NuclearSE_withStatus', 
                    header = T) %>% dplyr::select(gene_ID,Nuclear,nucStatus)
readChromDF <- fread('NlaIII_GM12878_output/readChromunities_wClosestGene_Batche4379.gz',nThread = 8)
colnames(readChromDF) <- c("ReadID","Cardinality","readLength","readQual","gene_ID",
                           "codingStatus","geneName","distToGene")

## Distribution of splicing efficiency in GM12878
ggplot(SE_df,aes(x = Nuclear)) + geom_histogram()

## less than 100KB away from gene
filteredReadChromDF <- readChromDF %>% filter(distToGene <= 1e5)

length(unique(readChromDF$ReadID))
length(unique(filteredReadChromDF$ReadID))
15984564*100/16099109 = 99.28

## Number of reads per cardinality
numReadsPerC <- as.data.frame(table(readChromDF %>% select(ReadID,Cardinality) %>% distinct() %>% .$Cardinality))
numReadsPerC <- numReadsPerC %>% mutate(perc = Freq*100/sum(Freq))
numReadsPerC$cumuPerc <- cumsum(numReadsPerC$perc)
numReadsPerC$Var1 <- as.numeric(as.character(numReadsPerC$Var1))
ggplot(numReadsPerC, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  xlim(0,30) + ggtitle("Read concatemers") +
  xlab("Cardinality")

ggplot(numReadsPerC, aes(x = Var1, y = cumuPerc)) +
  geom_point() + geom_line() +
  xlim(0,30) + ggtitle("Percent of data") +
  xlab("Cardinality") + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))

## Read length distribution
ggplot(readChromDF %>% filter(Cardinality <= 20), 
       aes(x = as.factor(Cardinality), y = readLength)) +
  geom_boxplot(outlier.shape = NA) +
  ggtitle("Read length distribution") +
  scale_y_continuous(trans = "log10")


## Simulating background/ null distribution of splicing efficiency
readCard_distr <- sample(2:30, size = 1e5, prob = numReadsPerC$perc[1:29]/100, replace = TRUE)
simReadCardStats <- as.data.frame(table(readCard_distr))
simReadCardStats$readCard_distr <- as.integer(as.character(simReadCardStats$readCard_distr))
ggplot(simReadCardStats, aes(x = readCard_distr, y = Freq)) +
  geom_bar(stat = "identity") +
  ggtitle("Simulated read concatemers") +
  xlab("Cardinality") +
  xlim(0,30)

## Randomly sample from genes to get a read representing c genes and get SE for those genes
sample_SE_forPoreC_reads <- function(cardin,rID){
  read_rID <- sample(SE_df$Nuclear,cardin)
  stats_rID <- as.vector(summary(read_rID))
  disp_rID <- max(stats_rID[6]-stats_rID[3],stats_rID[3]-stats_rID[1])
  mean_rID <- stats_rID[4]
  median_rID <- stats_rID[3]
  df <- data.frame(readID = rID, cardinality = cardin, 
             meanSE = mean_rID,medianSE = median_rID,
             dispersion = disp_rID)
  return(df)
}

set.seed(10)
null_DF <- do.call("rbind",lapply(1:20, 
                        function(i) do.call("rbind", lapply(1:simReadCardStats$Freq[i], 
                        function(rID) sample_SE_forPoreC_reads(rID, cardin = simReadCardStats$readCard_distr[i]) ))))


null_DF$cardinality <- as.factor(as.numeric(null_DF$cardinality))
ggplot(null_DF,aes(x = cardinality, y = meanSE)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.7,1)

ggplot(null_DF,aes(x = cardinality, y = medianSE)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.7,1)

ggplot(null_DF,aes(x = cardinality, y = dispersion)) +
  geom_boxplot(outlier.shape = NA)

## same number of reads per cardinality values
set.seed(10)
null_sameDF <- do.call("rbind",lapply(1:20, 
                                  function(i) do.call("rbind", lapply(1:500, 
                                                                      function(rID) sample_SE_forPoreC_reads(rID, cardin = simReadCardStats$readCard_distr[i]) ))))


null_sameDF$cardinality <- as.factor(as.numeric(null_sameDF$cardinality))
ggplot(null_sameDF,aes(x = cardinality, y = meanSE)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.7,1)

ggplot(null_sameDF,aes(x = cardinality, y = medianSE)) +
  geom_boxplot(outlier.shape = NA) + ylim(0.7,1)

ggplot(null_sameDF,aes(x = cardinality, y = dispersion)) +
  geom_boxplot(outlier.shape = NA)

## Real splicing efficiency values for Pore-C data
## Correlate with splicing -------
## very basic QC 
filteredReadChromDF_wSE <- left_join(sample_n(filteredReadChromDF, 1e6),SE_df,by = "gene_ID")
ggplot(filteredReadChromDF_wSE,aes(x = Nuclear)) + geom_histogram()

## Unbiased by gene expression (rather counts)
ggplot(filteredReadChromDF_wSE %>% select(geneName, Nuclear) %>% distinct(),
       aes(x = Nuclear)) + geom_histogram()

## get mean SE
combinedDF <- filteredReadChromDF_wSE %>% 
  group_by(ReadID) %>% mutate(meanSE=mean(Nuclear,na.rm = T), 
                              missValues = sum(is.na(Nuclear))/n())

combined_mod <- combinedDF %>% dplyr::select(ReadID,Cardinality,codingStatus,readLength,readQual,
                                             meanSE,missValues) %>% distinct()

## Filter low cardinality and fewer missing values for splicing efficiency:
combined_mod_filt <- combined_mod %>% filter(!is.na(meanSE) & 
                                               Cardinality <= 20 & missValues <= 0.5)
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

realDataCardDist <- as.data.frame(table(combined_mod_filt$Cardinality))
ggplot(realDataCardDist, aes(x = Var1,y = Freq)) +
  geom_bar(stat = "identity") + xlab("Cardinality")


## checking to see if missing values favor one cardinality over another
## This may be interesting to see which + why these genes are missing
readsPerMissingValueBin <- combinedDF %>% ungroup() %>% select(Cardinality,missValues) %>% 
  mutate(missBin = cut(missValues,breaks = seq(0,1,0.25),include.lowest = TRUE)) %>% 
  group_by(Cardinality,missBin) %>% summarize(numReads = n(),.groups = "keep") %>%
  ungroup() %>% group_by(Cardinality) %>% mutate(Perc = numReads/sum(numReads)) %>%
  as.data.frame()

ggplot(readsPerMissingValueBin %>% filter(Cardinality <= 9),
       aes(x = missBin, y = Perc, fill = as.factor(Cardinality))) +
    geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none") + 
  xlab("Missing SE values (%)") 




