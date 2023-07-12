#!/bin/R
# R version 4.1.2
# Synergy formulation of Pore-C data

# Setup -----
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

# Functions from new chromunity -----

csv2gr = function(path = NULL, col_names = NULL, save_path = NULL, prefix = "NlaIII_this_sample", mc.cores = 5, verbose = TRUE){
  
  if(is.null(path)){
    stop("Need a valid path to all Pore-C csvs.")
  }
  
  all.paths = data.table(file_path = dir(path, all.files = TRUE, recursive = TRUE, full = TRUE))[grepl("*fragment_alignments.csv.gz*", file_path)]
  
  if(nrow(all.paths) == 0){
    stop("No valid files files with suffix pore_c.csv found.")
  } 
  
  if(verbose){"Beginning to read csv files"}
  
  if(is.null(col_names)){col_names = c("read_name", "chrom", "start", "end", "pass_filter")}
  
  csv.list = pbmclapply(1:nrow(all.paths), function(k){
    csv.al = fread(all.paths[k]$file_path)[, .(read_name, chrom, start, end, pass_filter)]
    csv.al = as.data.table(csv.al)
    csv.al = csv.al[pass_filter ==  TRUE]  
    return(csv.al)
  }, mc.cores = mc.cores)
  
  csv.dt = rbindlist(csv.list, fill = TRUE)
  
  gc()
  
  csv.dt[, read_idx := .GRP, by = read_name]
  uni.dt = unique(csv.dt[, .(read_name, read_idx)])
  uni.dt[, fr.read_name := .N, by = read_name]
  uni.dt[, fr.read_idx := .N, by = read_idx]
  if (!any(uni.dt$fr.read_name > 1)){
    if(!any(uni.dt$fr.read_idx > 1)){
      rm(uni.dt)
    } else {
      message("Duplicate read_idx found, check this field in the output")
    }
  } else {
    message("Duplicate read_name found, make all files belong to a single, unique sample")
  }
  
  csv.gr = dt2gr(csv.dt)
  
  if (!is.null(save_path)){
    saveRDS(csv.gr, paste0(save_path, "/", prefix, ".rds"))
  }
  return(csv.gr)
}

# Load libraries -------
library(data.table)
library(dplyr)
library(rtracklayer)
library(skitools)
library(chromunity)
library(MASS)
library(gUtils)
library(parallel)
library(pbmcapply)
library(arrow)


# Import data and divide into training and testing
## Generating gRanges file from CSV
testFile <- csv2gr('../v0_poreC_explore/',mc.cores = 2)
toRemove <- c(grep("_",seqlevelsInUse(testFile),value = T),"chrEBV","chrM")
testFile <- dropSeqlevels(testFile,toRemove,pruning.mode = "coarse")
testFile$cid <- testFile$read_idx

## Customizing input
testFile <- fread('../v0_poreC_explore/fragFile_sorted_wClosestGene.bed.gz')[,c(1:4,11)]
colnames(testFile) <- c("seqnames","start","end","Read","Gene")
testFile <- testFile %>% tidyr::separate(Read,c("R","read_idx","C","card"),sep = ":|_") %>% select(-c(R,C))
testFile <- testFile %>% mutate_at("read_idx",as.integer)
testFile <- dt2gr(testFile)

## Isolate a chromosome or divide gRanges into multiple chromosomes
# test_chr18 <- testFile %Q% (seqnames == "chr18")
# 
# all_concatemers = unique(test_chr18$read_idx)
# set.seed(125)
# training_concatemers = sample(all_concatemers, length(all_concatemers)/2)
# testing_concatemers = setdiff(all_concatemers, training_concatemers)
# 
# test_chr18_training = test_chr18 %Q% (read_idx %in% training_concatemers)
# test_chr18_testing = test_chr18 %Q% (read_idx %in% testing_concatemers)
# 
# ##
# test_chr18_training$cid = test_chr18_training$read_idx
# test_chr18_testing$cid = test_chr18_testing$read_idx
# 
# ##
# test_chr18_sliding_chrom = chromunity(concatemers = test_chr18_training,piecewise = F,
#                                       resolution = 1e4, window.size = 2e6, mc.cores = 1)

debug(sliding_window_chromunity)
this_sliding_chrom = sliding_window_chromunity(concatemers = testFile, resolution = 5e5, window.size = 2e7, 
                                               take_sub_sample = FALSE, chr  = "chr8", subsample.frac = NULL, mc.cores = 1)
# this_sliding_chrom = sliding_window_chromunity(concatemers = testFile, resolution = 5e4, window.size = 2e6, 
#                                                take_sub_sample = TRUE, chr  = "chr8", subsample.frac = 0.5, mc.cores = 1)

## 

## Generate covariates ------
## Frag counts
#frags <- read_parquet("../v0_poreC_explore/DpnII_GRCh38.vd.fragments.parquet") %>% dt2gr()
frags = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/Nlaiii.frags.hg38.rds"))) %>% dt2gr()

## GC content
gc5b = readRDS(gzcon(file("https://mskilab.s3.amazonaws.com/chromunity/gc.38.rds")))
cov_list = list(gc5b, frags)

## Specify what kind of covariate it is. Score will be aggregated over bins while the number of intervals will be calculated otherwise.
names(cov_list) <- c("score:gc.cov", "interval:frag.cov")
gc_frag_cov = covariate(name = c("gc", "frag"), type = c("numeric", "interval"), field = c("score", NA), data = cov_list)


## Annotate
## concatemers used for testing
this.chrom = gr2dt(this_sliding_chrom$concatemers)

## Concatemers for testing
this_gr_testing = dt2gr(gr2dt(testFile)[!read_idx %in% unique(this.chrom$read_idx)]) 

debug(annotate)
annotated_chrom = annotate(binsets = this_sliding_chrom$binsets,
                           k = 5,
                           concatemers = this_gr_testing,
                           covariates = gc_frag_cov, resolution = 5e5,
                           mc.cores = 3) 

