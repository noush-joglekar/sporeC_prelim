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
testFile <- fread('NlaIII_GM12878_output_byChr/NlaIII_GM12878_chr8.gz')[,c(1:3,13,15)]
colnames(testFile) <- c("seqnames","start","end","Gene","read_idx")
testFile <- testFile %>% mutate_at("read_idx",as.integer)
testFile <- dt2gr(testFile)


## Concatemer shuffling function (from Jameson):
shuffle_concatemers = function(concatemers, contact_matrix) {
  A = contact_matrix$mat %>% as.matrix
  rownames(A) <- NULL
  colnames(A) <- NULL
  A[cbind(1:nrow(A), 1:nrow(A))] = 0
  A = A + t(A)
  An = A
  An = round(1+10*An/min(An[An>0]))  ##can you remove this background 1 value? does this change anything. peculiar     An = round(1+10*An/min(An[An>0]))
  edges = as.data.table(which(An != 0, arr.ind = TRUE))[ , val := An[which(An!=0)]][, .(row = rep(row, val), col = rep(col, val))]
  ##
  G = graph.edgelist(edges[, cbind(row, col)])
  
  concats.dt = gr2dt(concatemers)
  concat.counts = concats.dt[, new.count := .N, by='read_idx']
  card = unique(concat.counts[, .(read_idx, new.count)])
  this.steps = sum(card$new.count)
  
  RW = random_walk(G, start = 200, steps = sum(card$new.count)) %>% as.numeric
  ##
  rm(G)
  rm(edges)
  ##gc()
  out = contact_matrix$gr[RW]%>% gr2dt()
  out$read_idx = card[, rep(read_idx, new.count)]
  ##out[, bid := all.bid[nr]]
  ##out[, cid := read_idx]
  return(out)
}


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




# RE chromunity tutorial ------

## create test file
## their file
# chr8_parq_gr = readRDS("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/tutorial_inputs/chr8.single.run.sub.rds")


## getting E-P annotations
chain19to38 = rtracklayer::import.chain("/gpfs/commons/groups/imielinski_lab/projects/PoreC/db/tutorial_inputs/hg19ToHg38.over.chain")
chmm_19 = import("https://mskilab.s3.amazonaws.com/chromunity/wgEncodeBroadHmmGm12878HMM.bed")
chmm = chmm_19 %>% rtracklayer::liftOver(chain19to38)%>% grl.unlist %Q% (!duplicated(grl.ix))
chmm_ep = dt2gr(gr2dt(chmm)[grepl("Promoter", name) | grepl("Enhancer", name)])

## exons:
gtf <- import('/gpfs/commons/groups/gursoy_lab/ajoglekar/Support/References/Human/GRCh38/GENCODE/gencode.v42.annotation.gtf')
exomeLines <- gtf %Q% (type == "exon")
targets_exome = gr.reduce(exomeLines+resolution)
### 

## resolution of the pad
resolution = 2.5e4 #4e4
tiles = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), resolution)

## Targets
this_EP = (chmm_ep %Q% (grepl("Active_Promoter|Strong_Enhancer", name)))
targets_EP = gr.reduce(this_EP+resolution)

## Subsample
chr8_parq_gr <- testFile ## ***
chr8_parq_gr = chr8_parq_gr %&% targets_exome #targets_EP

## Subsampling concatemers
set.seed(125)
all_concatemers = unique(chr8_parq_gr$read_idx)
training_concatemers = sample(all_concatemers, length(all_concatemers)/2)
testing_concatemers = setdiff(all_concatemers, training_concatemers)

this_gr_training = chr8_parq_gr %Q% (read_idx %in% training_concatemers)
this_gr_testing = chr8_parq_gr %Q% (read_idx %in% testing_concatemers)

##
this_gr_training$cid = this_gr_training$read_idx
this_gr_testing$cid = this_gr_testing$read_idx

## Running Chromunity
this_re_chrom = re_chromunity(concatemers = this_gr_training, windows = targets_exome, #targets_EP
                              piecewise = FALSE, shave = TRUE, resolution = 2.5e4, mc.cores = 5) #1e4 4e4


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

gc_frag_cov

### Bin annotation
set.seed(198)
annotated_re_chrom = chromunity::annotate(binsets = this_re_chrom$binsets,
                              k = 3,
                              concatemers = this_gr_testing,
                              covariates = gc_frag_cov, resolution = 4e4, #2.5e4
                              mc.cores = 5) 
head(annotated_re_chrom)


## background sets
set.seed(198)
back_re_gr = gr2dt(dt2gr(re_background(binsets = this_re_chrom$binsets, n = 1000,
                                       resolution = 2.5e4)))

upper.bound = as.data.table(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), keep.rownames = T) 
setkeyv(back_re_gr, c("seqnames", "start"))
back_re_gr = back_re_gr[!bid %in% back_re_gr[width < (2.5e4-1)]$bid]
back_re_gr = gr2dt(gr.reduce(dt2gr(back_re_gr), by = "bid"))
back_re_gr$bid <- as.factor(back_re_gr$bid)
back_re_gr = merge(back_re_gr, upper.bound, by.x = "seqnames", by.y = "V1", all.x = T, allow.cartesian = T)
back_re_gr = back_re_gr[end < V2][start < V2]
back_re_gr[, overall.cardinality := .N, by = bid]
back_re_gr = back_re_gr[overall.cardinality > 1]
back_re_gr = dt2gr(back_re_gr)
head(back_re_gr)


set.seed(198)
annotated_re_back = chromunity::annotate(binsets = back_re_gr,
                             k = 3,
                             concatemers = this_gr_testing,
                             covariates = gc_frag_cov, resolution = 2.5e4,
                             mc.cores = 5) 

## Removing edge cases with no counts
annotated_re_back = annotated_re_back[!bid %in% annotated_re_back[, .(sum(count)), by = bid][V1 == 0]$bid]
head(annotated_re_back)



## Train model
set.seed(198)
back_re_model = fit(annotated_re_back)

annotated_re_chrom = sscore(annotated_re_chrom, model = back_re_model) 
head(annotated_re_chrom)


## Synergy model

set.seed(198)
this_synergy = na.omit(synergy(binsets = this_re_chrom$binsets,
                               annotated.binsets = annotated_re_chrom, model = back_re_model))


this_synergy$fdr = signif(p.adjust(this_synergy$p, "BH"), 2)

reads_wBIds <- chr8_parq_gr %*% this_re_chrom$binsets

sigReads <- reads_wBIds %Q% (bid %in% unique(this_synergy %>% filter(fdr <= 0.05) %>% .$bid) )
sigSyns <- gr2dt(sigReads) %>% dplyr::select(Gene,bid) %>% distinct()


nonSigReads <- reads_wBIds %Q% (!(bid %in% this_synergy$bid))
nonSigBids <- gr2dt(nonSigReads) %>% dplyr::select(Gene,bid) %>% distinct()


## Get shuffled concatemers
inMat <- read.csv('tmpDir/projMat_card3_chrom8.csv')
shuffled <- shuffle_concatemers()

library(ggVennDiagram)

geneNames <- list("SigGenes" = sigSyns$Gene, "NonSigGenes" = nonSigBids$Gene)
ggVennDiagram(geneNames, category.names = c("Sig","NonSig")) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "bottom")


binNames <- list("SigBins" = sigSyns$bid, "NonSigBins" = nonSigBids$bid)
ggVennDiagram(binNames, category.names = c("Sig","NonSig")) + 
  scale_fill_gradient(low = "#fadede", high = "#db4848") +
  theme(legend.position = "bottom")

## Bring in splicing efficiency

SE_genes <- read.table('Gene_NucSE',sep = "\t")
colnames(SE_genes) <- c("GeneID","Gene","SE")

sigSyns_wSE <- left_join(sigSyns,SE_genes) %>% filter(!(Gene %in% nonSigBids$Gene)) %>% mutate(Status = "Sig")
nonSig_wSE <- left_join(nonSigBids,SE_genes)  %>% filter(!(Gene %in% sigSyns$Gene)) %>% mutate(Status = "NonSig")

allSE <- rbind(sigSyns_wSE,nonSig_wSE)

allSE$bid <- as.factor(allSE$bid)
ggplot(allSE, aes(x = Status, y = SE, fill = Status)) +
  geom_boxplot() +ylim(0.75,1.03) + 
  ggsignif::geom_signif(comparisons = list(c("Sig","NonSig")),
                        y_position = 0.98,map_signif_level = T,tip_length = 0.01)
