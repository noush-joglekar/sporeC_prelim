#!/bin/R
# R version 3.6.0
# Synergy formulation of Pore-C data: Tutorial

scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

library(data.table)
library(dplyr)
library(rtracklayer)
library(skitools)
library(chromunity)
library(MASS)
library(chromunity)
library(gUtils)
library(parallel)


## following tutorial
this_gr = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/chr8.single.run.sub.rds")))
this_gr$cid = this_gr$read_idx


this_sliding_chrom = sliding_window_chromunity(concatemers = this_gr, chr = "chr8",
                                               take_sub_sample = TRUE,subsample.frac = 0.5,
                                               resolution = 5e4, window.size = 2e6, mc.cores = 1)


## Frag counts
## frags = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/tutorial_inputs/Nlaiii.frags.hg38.rds") %>% dt2gr()
frags = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/Nlaiii.frags.hg38.rds"))) %>% dt2gr()

## GC content
## gc5b = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/tutorial_inputs/gc.38.rds")
gc5b = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/gc.38.rds")))

## 
cov_list = list(gc5b, frags)

## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
names(cov_list) <- c("score:gc.cov", "interval:frag.cov")


gc_frag_cov = covariate(name = c("gc", "frag"), type = c("numeric", "interval"), field = c("score", NA), data = cov_list)


## concatemers used for testing
this.chrom = gr2dt(this_sliding_chrom$concatemers)

## Concatemers for testing
this_gr_testing = dt2gr(gr2dt(this_gr)[!read_idx %in% unique(this.chrom$read_idx)]) 

annotated_chrom = annotate(binsets = this_sliding_chrom$binsets,
                           k = 5,
                           concatemers = this_gr_testing,
                           covariates = gc_frag_cov, resolution = 5e4,
                           mc.cores = 3) 


## Background set on called chromunities
set.seed(198)
back_gr = gr2dt(dt2gr(background(binsets = this_sliding_chrom$binsets, n = 1000,
                                    resolution = 5e4)))

## Getting seqlengths of each chromosome
upper.bound = as.data.table(hg_seqlengths(genome = "../v0_poreC_explore/hg38.chromSizes"), keep.rownames = T) 
setkeyv(back_gr, c("seqnames", "start"))
back_gr = back_gr[!bid %in% back_gr[width < (5e4-1)]$bid]
back_gr = gr2dt(gr.reduce(dt2gr(back_gr), by = "bid"))
back_gr$bid <- as.factor(back_gr$bid)
back_gr = merge(back_gr, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
back_gr = back_gr[end < V2][start < V2]
back_gr[, overall.cardinality := .N, by = bid]
back_gr = back_gr[overall.cardinality > 1]
back_gr = dt2gr(back_gr)

## Annotate the background bins
set.seed(198)
annotated_back = annotate(binsets = back_gr,
                          k = 5,
                          concatemers = this_gr_testing,
                          covariates = gc_frag_cov, resolution = 5e4,
                          mc.cores = 3) 

## Removing edge cases with no counts
annotated_back = annotated_back[!bid %in% annotated_back[, .(sum(count)), by = bid][V1 == 0]$bid]

## Training the first model 
set.seed(198)
back_model = fit(annotated_back)

annotated_chrom = sscore(annotated_chrom, model = back_model) 
head(annotated_chrom)

## Second model to test for synergy
set.seed(198)
this_synergy = na.omit(synergy(binsets = this_sliding_chrom$binsets,
                               annotated.binsets = annotated_chrom, model = back_model))

this_synergy$fdr = signif(p.adjust(this_synergy$p, "BH"), 2)




