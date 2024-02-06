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

hcc <- read.table('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/v1.evaluateExpectedVersusInteresting_NlaIII_HCC1954/dfs_chr8/diffContactReads_chr8.tab.gz',
header = TRUE)
hcc <- as.data.frame(hcc) %>% select(chr,binStart,binEnd,binID,geneName,readID)

colnames(hg) <- colnames(hcc) <- c("seqnames","start","end","binID","Gene","readID")

hg_gr <- dt2gr(hg)
hcc_gr <- dt2gr(hcc)
view_range = GRanges('chr8:1000001-141000000')


plot(c(gTrack(split(hg_gr, hg_gr$readID) %>% unname, height = 10, name = 'HG002'), 
gTrack(split(hcc_gr, hg_gr$readID) %>% unname, height = 10, name = 'HCC1954')), view_range)
