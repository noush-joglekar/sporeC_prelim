#!/bin/R
# R version 4.1.2
# Synergy formulation of Pore-C data

## Setup -----
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/'
setwd(workingDir)

## Load libraries -------
library(chromunity)
library(dplyr,quietly = T,warn.conflicts = F)
library(tidyr,quietly = T,warn.conflicts = F)
library(tibble,quietly = T,warn.conflicts = F)
library(rtracklayer,quietly = T,warn.conflicts = F)
library(skitools,quietly = T,warn.conflicts = F)
library(MASS,quietly = T,warn.conflicts = F)
library(pbmcapply,quietly = T,warn.conflicts = F)
library(arrow,quietly = T,warn.conflicts = F)

print("All libraries loaded, reading in input now")

## Input ------
# Import data and divide into training and testing
## Customizing input
args <- commandArgs(trailingOnly = TRUE)
#args <- c('NlaIII_GM12878_output_byChr/NlaIII_GM12878_chr2.gz','chr2','synergyRun_byChrom/')
inFragFile <- args[1]
chrID <- args[2]
outputDir <- args[3]

print(paste0("Processing promethION file for ",chrID))

inputFile <- fread(inFragFile,nThread=3)[,c(1:3,13,15)]
colnames(inputFile) <- c("seqnames","start","end","Gene","read_idx")
inputFile <- inputFile %>% mutate_at("read_idx",as.integer)
inputFile <- dt2gr(inputFile)

inputFile$cid = inputFile$read_idx

print("Getting sliding window chromunities")
if(!file.exists(paste0(outputDir,"slidingWinChrom_",chrID,".rds"))){
  this_sliding_chrom = sliding_window_chromunity(concatemers = inputFile, resolution = 2e4, ## resolution = 5e4
                                                 window.size = 2e6, take_sub_sample = TRUE, 
                                                 chr  = chrID, subsample.frac = 0.5, mc.cores = 1)
  saveRDS(this_sliding_chrom,file = paste0(outputDir,"slidingWinChrom_",chrID,".rds"))
} else {
  print("Sliding window gRanges file exists -- moving on")
  this_sliding_chrom <- readRDS(paste0(outputDir,"slidingWinChrom_",chrID,".rds"))
}


## Covariates ------
print("Encoding covariates: NlaIII frags and GC content")

## GC content
if(!file.exists(file.path(outputDir,'/GC_fragCov.Rds'))){
  ## frags
  frags = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/Nlaiii.frags.hg38.rds"))) %>% dt2gr()
  
  gc5b = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/gc.38.rds")))
  cov_list = list(gc5b, frags)

  ## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
  gc_frag_cov = covariate(name = c("gc", "frag"), type = c("numeric", "interval"), field = c("score", NA), data = cov_list)
  
  saveRDS(gc_frag_cov,file = file.path(outputDir,'/GC_fragCov.Rds'))
} else {
  print("GC content + RE fragment covariate file -- moving on")
  gc_frag_cov <- readRDS(file.path(outputDir,'/GC_fragCov.Rds'))
}

print("Making GRanges of test set")
this.chrom = gr2dt(this_sliding_chrom$concatemers)
this_gr_testing = dt2gr(gr2dt(inputFile)[!read_idx %in% unique(this.chrom$read_idx)])

print("Annotating the sliding window chromunity")
if(!file.exists(paste0(outputDir,"annotatedChrom_",chrID,".rds"))){
  annotated_chrom = chromunity::annotate(binsets = this_sliding_chrom$binsets,
                                         k = 5,
                                         concatemers = this_gr_testing,
                                         covariates = gc_frag_cov, resolution = 2e4,
                                         mc.cores = 1)
  saveRDS(annotated_chrom,paste0(outputDir,"annotatedChrom_",chrID,".rds"))
} else {
  print("Annotated test set gRanges file exists -- moving on")
  annotated_chrom <- readRDS(paste0(outputDir,"annotatedChrom_",chrID,".rds"))
}

print("Building the background set")
## Background ------
set.seed(198)
if(!file.exists(paste0(outputDir,"slidingWindow_background_",chrID,".rds"))){
  back_gr = sliding_window_background(chromosome= chrID, binsets = this_sliding_chrom$binsets,
                                      resolution = 2e4)
  back_gr[, V1 := NULL]
  back_gr = na.omit(back_gr)
  ## Adding few filters to remove outlier simulation, can be customized
  ## Removing bins less than resolution and lying out of bounds
  
  ## Getting seqlengths of each chromosome
  upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T) 
  setkeyv(back_gr, c("seqnames", "start"))
  back_gr = back_gr[!bid %in% back_gr[width < (2e4-1)]$bid]
  back_gr = gr2dt(gr.reduce(dt2gr(back_gr), by = "bid"))
  back_gr$bid <- as.factor(back_gr$bid)
  back_gr = merge(back_gr, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
  back_gr = back_gr[end < V2][start < V2]
  back_gr[, overall.cardinality := .N, by = bid]
  back_gr = back_gr[overall.cardinality > 1]
  back_gr = dt2gr(back_gr)
  saveRDS(back_gr,paste0(outputDir,"slidingWindow_background_",chrID,".rds"))
} else {
  print("Background sliding window file exists -- moving on")
  back_gr <- readRDS(paste0(outputDir,"annotatedChrom_",chrID,".rds"))
}

print("Annotating the background set")
set.seed(198)

if(!file.exists(paste0(outputDir,"annotatedBG_",chrID,".rds"))){
  annotated_back = chromunity::annotate(binsets = back_gr,
                                        k = 5,
                                        concatemers = this_gr_testing,
                                        covariates = gc_frag_cov, resolution = 2e4,
                                        mc.cores = 1) 
  annotated_back = annotated_back[!bid %in% annotated_back[, .(sum(count)), by = bid][V1 == 0]$bid]
  saveRDS(annotated_back,paste0(outputDir,"annotatedBG_",chrID,".rds"))
} else {
  print("Annotated bg gRanges file exists -- moving on")
  annotated_back <- readRDS(paste0(outputDir,"annotatedBG_",chrID,".rds"))
}

print("Fit 1: Background model")
## Fit #1 Background ------
set.seed(198)
back_model = chromunity::fit(annotated_back)

annotated_chrom = sscore(annotated_chrom, model = back_model)

print("Fit 2: Test model")
## Fit #2 Obtain synergies ------
set.seed(198)
this_synergy = na.omit(synergy(binsets = this_sliding_chrom$binsets,
                               annotated.binsets = annotated_chrom, model = back_model))

this_synergy$fdr = signif(p.adjust(this_synergy$p, "BH"), 2)

print("Obtaining significant synergies")

## Sig synergy: -----
sigSyn <- this_synergy %>% filter(fdr <= 0.1)
print(paste("There are",nrow(sigSyn),"synergistic chromunities for",chrID))


print("Saving output")
write.table(gr2dt(this_synergy), paste0(outputDir,"allSynergies_",chrID,".tab"),
            sep = "\t",quote = F, row.names = F, col.names = T)
print("Done!")
