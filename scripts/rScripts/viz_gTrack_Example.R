library(chromunity)
library(tidyr)
library(dplyr)
library(gUtils)
library(gTrack)


simplicies = readRDS('/gpfs/commons/groups/imielinski_lab/home/jorvis/Projects/testing/db/random_simplex_snippet/simplicial_snippet.rds')
view_range = GRanges('chr8:126000000-129000000')
plot(gTrack(split(simplicies, simplicies$cid) %>% unname, height = 10, name = 'simplicies'), view_range)


## Cancer cell line data
hg <- read.table('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/v1.evaluateExpectedVersusInteresting_NlaIII_HG002/dfs_chr8/diffContactReads_chr8.tab.gz',
header = TRUE)
hg <- as.data.frame(hg) %>% select(chr,binStart,binEnd,binID,geneName,readID)

# hcc <- read.table('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/v1.evaluateExpectedVersusInteresting_NlaIII_HCC1954/dfs_chr8/diffContactReads_chr8.tab.gz',
# header = TRUE)
hcc <- read.table('diffAnalysis/diffContactReads_chr11_HCC1954_compGM.tab.gz',header = T)
hcc <- as.data.frame(hcc) %>% select(chr,binStart,binEnd,binID,geneName,readID)

colnames(hg) <- colnames(hcc) <- c("seqnames","start","end","binID","Gene","readID")

hg_gr <- dt2gr(hg)
hcc_gr <- dt2gr(hcc)
view_range = GRanges('chr8:1000001-141000000')

plot(c(gTrack(split(hg_gr, hg_gr$readID) %>% unname, height = 10, name = 'HG002'), 
gTrack(split(hcc_gr, hg_gr$readID) %>% unname, height = 10, name = 'HCC1954')), view_range)

## GM12878 data
# gm <- read.table('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/v1.evaluateExpectedVersusInteresting_NlaIII_GM12878_2/dfs_chr8/diffContactReads_chr8.tab.gz',
#                  header = TRUE)
gm <- read.table('diffAnalysis/diffContactReads_chr11_GM12878_compHCC.tab.gz', header=TRUE)
gm <- as.data.frame(gm) %>% select(chr,binStart,binEnd,binID,geneName,readID)
colnames(gm) <- c("seqnames","start","end","binID","Gene","readID")
gm_gr <- dt2gr(gm)
view_range = GRanges('chr11:500001-135000000')


plot(gTrack(split(gm_gr, gm_gr$readID) %>% unname, height = 10, name = 'GM12878'), view_range)

## Random stats 
hg_genes = unique(hg$Gene)
hcc_genes = unique(hcc$Gene)
both_genes = intersect(hg_genes,hcc_genes)
