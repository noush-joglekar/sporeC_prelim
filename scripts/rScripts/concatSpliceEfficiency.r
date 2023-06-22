#!/bin/R

## Multiway interaction project
## Concatenate splicing efficiency data across replicates for GM12878 various cellular fractions
## See if we can divide genes into low, medium, high SE

scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2023_03_06_GM12878_cellularFractionData/v0.evaluateSplicingEfficiency/'
setwd(workingDir)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(enrichR)


samples <- c("WholeCell","Nuclear","Cytosolic")
replicates <- c("Rep1","Rep2")

geneInfo <- read.table('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_explore/geneCoordsAndNames_sorted.bed',
                       sep="\t", header = F)[,c(5:7)]
colnames(geneInfo) <- c("gene_ID","bioType","gene_name")
geneInfo <- geneInfo %>% mutate_all(~gsub(";","",.))

sampleIDs <- sapply(replicates, function(r) paste(samples,r, sep = "_"))

seDF <- NULL

for(s in samples){
  seFiles <- list.files(paste0('../GM12878_',grep(s,sampleIDs, value = T),'/spliceq/'), full.names = T)
  for (f in seFiles){
    sampleID = paste(unlist(strsplit(basename(f),"_|[.]"))[2:3],collapse = "_")
    contents <- fread(f) %>% select(chr,gene_ID,transcript_ID,intron_ID, score)
    contents$SampleID = sampleID
    seDF[[sampleID]] <- contents
  }
}

seBound <- do.call('rbind',seDF)
seBound <- left_join(seBound,geneInfo,by = "gene_ID")

seBound_avgOverRep <- seBound %>% 
  separate(SampleID, into = c("Sample","Rep"),sep = "_") %>%
  group_by(gene_ID,gene_name,intron_ID,Sample,bioType) %>%
  summarise(rmS = mean(score,na.rm = T)) %>% ungroup() %>% 
  group_by(gene_ID, Sample) %>% mutate(imS = mean(rmS,na.rm = T), imSD = sd(rmS))

## avg per replicate for nuclear
ggplot(seBound_avgOverRep %>% filter(Sample == "Nuclear" & bioType %in% c("lncRNA","protein_coding")),
       aes(x = gene_name, y = rmS, col = bioType)) +
  geom_point(alpha = 0.6, size = 0.1) +
  theme(axis.text.x = element_blank())

## avg per replicate and gene for nuclear
ggplot(seBound_avgOverRep %>% filter(Sample == "Nuclear" & bioType %in% c("lncRNA","protein_coding")),
       aes(x = gene_name, y = imS, col = bioType)) +
  geom_point(alpha = 0.6, size = 0.1) +
  theme(axis.text.x = element_blank()) 

## avg per replicate and gene for all as a sanity check. Nuclear always less than cytosolic and whole cell
ggplot(seBound_avgOverRep %>% filter(bioType %in% c("lncRNA","protein_coding")),
       aes(x = gene_name, y = imS, col = Sample)) +
  geom_point(alpha = 0.6, size = 0.1) +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~bioType)

lowSE_genes <- seBound_avgOverRep %>% filter(Sample == "Nuclear" & bioType == "protein_coding") %>% 
       filter(imS <= 0.25)

ggplot(lowSE_genes,
       aes(x = gene_name, y = rmS)) +
  geom_point(alpha = 0.6) + geom_line() +
  theme(axis.text.x = element_blank()) 

dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
if (websiteLive) {
  enriched <- enrichr(unique(lowSE_genes$gene_name), dbs)
}

enrichR::plotEnrich(enriched$GO_Biological_Process_2018)

## okay back to averaging over all introns of a gene, no harm done
seBound_avgPerGene <- seBound %>% 
  separate(SampleID, into = c("Sample","Rep"),sep = "_") %>%
  group_by(gene_ID,Sample) %>%
  summarise(mS = mean(score,na.rm = T))

seBound_avgPerGene_wVar <- seBound_avgPerGene %>% group_by(gene_ID) %>%
  mutate(v = var(mS)) %>%
  pivot_wider(names_from = Sample, values_from = mS)

hist(seBound_avgPerGene_wVar$v,breaks = 50) ## very little discordance

smoothScatter(seBound_avgPerGene_wVar$Nuclear,seBound_avgPerGene_wVar$Cytosolic,
              xlab = "Nuclear", ylab = "Cytosolic", main = "SE correlation")
smoothScatter(seBound_avgPerGene_wVar$WholeCell,seBound_avgPerGene_wVar$Cytosolic,
              xlab = "Cytosolic", ylab = "WholeCell", main = "SE correlation")

hist(seBound_avgPerGene_wVar$Nuclear,main = "Nuclear fraction", 
     xlab = "Splicing efficiency GM12878")

seBound_avgPerGene_wVar <- seBound_avgPerGene_wVar %>%
  filter(!is.na(Nuclear)) %>%
  mutate(nucStatus = case_when(Nuclear <= 0.25 ~ "Low",
                               Nuclear > 0.25 & Nuclear <= 0.75 ~ "Medium",
                               Nuclear > 0.75 ~ "High")) %>%
  rename(SE_var = v)

seBound_avgPerGene_wVar$nucStatus <- factor(seBound_avgPerGene_wVar$nucStatus, 
                                            levels = c("Low","Medium","High"))

ggplot(seBound_avgPerGene_wVar, 
       aes(x = nucStatus, y = Nuclear, fill = nucStatus)) +
  geom_boxplot() + theme_classic() 


write.table(seBound_avgPerGene_wVar, file = 'SE_AvgOverReps_NuclearSE_withStatus',
            sep = "\t",quote = F, row.names = FALSE, col.names = TRUE)


