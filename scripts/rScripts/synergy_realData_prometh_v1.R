#!/bin/R
# R version 4.1.2
# Synergy formulation of Pore-C data

# Setup -----
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/'
setwd(workingDir)

# Load libraries -------
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(rtracklayer)
library(skitools)
library(chromunity)
library(MASS)
library(gUtils)
library(parallel)
library(pbmcapply)
library(arrow)
library(edgeR)

# Import data and divide into training and testing
## Customizing input
testFile <- fread('NlaIII_GM12878_output_byChr/NlaIII_GM12878_chr19.gz')[,c(1:3,13,15)]
colnames(testFile) <- c("seqnames","start","end","Gene","read_idx")
testFile <- testFile %>% mutate_at("read_idx",as.integer)
testFile <- dt2gr(testFile)


testFile$cid = testFile$read_idx
this_sliding_chrom = sliding_window_chromunity(concatemers = testFile, resolution = 5e4, 
                                               window.size = 2e6, take_sub_sample = TRUE, 
                                               chr  = "chr19", subsample.frac = 0.5, mc.cores = 4)

## frags
frags = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/Nlaiii.frags.hg38.rds"))) %>% dt2gr()

## GC content
gc5b = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/gc.38.rds")))
cov_list = list(gc5b, frags)

## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
names(cov_list) <- c("score:gc.cov", "interval:frag.cov")
gc_frag_cov = covariate(name = c("gc", "frag"), type = c("numeric", "interval"), field = c("score", NA), data = cov_list)

this.chrom = gr2dt(this_sliding_chrom$concatemers)
this_gr_testing = dt2gr(gr2dt(testFile)[!read_idx %in% unique(this.chrom$read_idx)])

annotated_chrom = annotate(binsets = this_sliding_chrom$binsets,
                           k = 5,
                           concatemers = this_gr_testing,
                           covariates = gc_frag_cov, resolution = 5e4,
                           mc.cores = 3) 

## Background
set.seed(198)
back_gr = sliding_window_background(chromosome= "chr19", binsets = this_sliding_chrom$binsets, n = 1000,
                                    resolution = 5e4)
## Generating distributions
## Generating GRanges
back_gr[, V1 := NULL]
back_gr = na.omit(back_gr)
## Adding few filters to remove outlier simulation, can be customized
## Removing bins less than resolution and lying out of bounds

## Getting seqlengths of each chromosome
upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T) 
setkeyv(back_gr, c("seqnames", "start"))
back_gr = back_gr[!bid %in% back_gr[width < (5e4-1)]$bid]
back_gr = gr2dt(gr.reduce(dt2gr(back_gr), by = "bid"))
back_gr$bid <- as.factor(back_gr$bid)
back_gr = merge(back_gr, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
back_gr = back_gr[end < V2][start < V2]
back_gr[, overall.cardinality := .N, by = bid]
back_gr = back_gr[overall.cardinality > 1]
back_gr = dt2gr(back_gr)

set.seed(198)
annotated_back = annotate(binsets = back_gr,
                          k = 5,
                          concatemers = this_gr_testing,
                          covariates = gc_frag_cov, resolution = 5e4,
                          mc.cores = 3) 

annotated_back = annotated_back[!bid %in% annotated_back[, .(sum(count)), by = bid][V1 == 0]$bid]

## Fit #1
set.seed(198)
back_model = chromunity::fit(annotated_back)

annotated_chrom = sscore(annotated_chrom, model = back_model)
head(annotated_chrom)

## Obtain synergies #2
set.seed(198)
this_synergy = na.omit(synergy(binsets = this_sliding_chrom$binsets,
                               annotated.binsets = annotated_chrom, model = back_model))
## Synergy 2022-08-17 20:25:44: Scoring binsets

this_synergy$fdr = signif(p.adjust(this_synergy$p, "BH"), 2)
head(this_synergy[order(p)])

## Sig synergy: 
sigSyn <- this_synergy %>% filter(fdr <= 0.1)
sigConcatemers <- gr2dt(this_sliding_chrom$concatemers %Q% (chid %in% sigSyn$bid)) %>%
  mutate(Status = "Sig")

nonSig_ctrl <- this_synergy %>% filter(p > 0.1) %>% slice_sample(n = nrow(sigSyn))
nonSigConcatemers <- gr2dt(this_sliding_chrom$concatemers %Q% (chid %in% nonSig_ctrl$bid)) %>%
  mutate(Status = "NonSig")

## Bring in splicing efficiency:
SE_genes <- read.table('Gene_NucSE',sep = "\t")
colnames(SE_genes) <- c("GeneID","Gene","SE")

sigConcat_SE <- left_join(rbind(sigConcatemers,nonSigConcatemers) %>% 
                            dplyr::select(Gene,read_idx,Status),SE_genes,by = "Gene")
ggplot(sigConcat_SE %>% drop_na(), 
       aes(x = Status, y = SE, fill = Status)) +
  geom_boxplot(outlier.alpha = 0.6) +
  ylim(c(0.75,1))

## Bring in gene expression 
rep1 <- read.table('../2023_03_06_GM12878_cellularFractionData/GM12878_Nuclear_Rep1/featureCountsOut/GM12878_Nuclear_Rep1.gene.counts.txt',
                   header = TRUE)
rep2 <- read.table('../2023_03_06_GM12878_cellularFractionData/GM12878_Nuclear_Rep2/featureCountsOut/GM12878_Nuclear_Rep2.gene.counts.txt',
                   header = TRUE)
gEx <- left_join(rep1,rep2) %>% column_to_rownames("Geneid")
colnames(gEx) <- c("Rep1","Rep2")
norm_gEx <- cpm(gEx)
smoothScatter(log10(norm_gEx[,1]),log(norm_gEx[,2]), xlab = "Log norm CPM Rep1", ylab = "Log norm CPM Rep2")

norm_gEx_df <- as.data.frame(norm_gEx) %>% rownames_to_column("GeneID")
sigConcat_SE_gEx <- left_join(sigConcat_SE,norm_gEx_df,by = "GeneID") %>% 
  rowwise() %>% mutate(logCPM = log10(mean(c(Rep1,Rep2))))


ggplot(sigConcat_SE_gEx %>% drop_na(), 
       aes(x = Status, y = logCPM, fill = Status)) +
  geom_boxplot(outlier.alpha = 0.6)
