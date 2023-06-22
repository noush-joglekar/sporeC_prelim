#!/bin/R
# By Anoushka Joglekar 04.2023
# V2
## Reads in a sparse matrix containing pairwise gene x gene interactions
## from reads spanning multiple genes. 
## I am starting with protein-coding genes only 
## It normalizes the full matrix using some version of expected versus observed

library(Matrix)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

# Reading input ------
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

## read in concatemers
readChromDF <- fread('../v0_poreC_explore/readChromunities_wClosestGene.gz')
colnames(readChromDF) <- c("ReadID","Cardinality","gene_ID","geneType","geneName","distToGene")

## read in gene coordinates
geneCoords <- read.table('../v0_poreC_explore/geneCoordsAndNames_sorted.bed')
colnames(geneCoords) <- c("chr","start","end","strand","geneID",
                          "s1","type","s2","geneName","s3")
geneCoords <- geneCoords %>% dplyr::select(-c(s1,s2,s3))

## read in chromosome sizes
chromSizes <- read.table('../v0_poreC_explore/hg38.chromSizes')
colnames(chromSizes) <- c("chr","size")

## Bin the genes into chromosome blocks 
binWidth = 1e7 #(10 MB)
bG <- list()
bG <- lapply(1:(nrow(chromSizes)-1), function(cZ) {
  geneCoords %>% filter(chr == chromSizes[cZ,1]) %>% 
  mutate_at(c("start","end"),as.integer) %>%
  mutate( coordStart = cut( start, breaks = seq(1,chromSizes[cZ,2],binWidth) ),
          coordEnd = cut( end, breaks = seq(1,chromSizes[cZ,2],binWidth) )) %>%
  group_by(coordStart) %>% mutate(binID = as.numeric(str_extract(coordStart,"(?<=,)[^\\]]+"))/binWidth)})
binnedGenes <- do.call('rbind',bG) %>% ungroup()

binnedGeneDF <- binnedGenes %>% select(chr,binID,geneID,geneName) %>% 
  group_by(chr,binID) %>% mutate(gID = cur_group_id()) %>%
  unite(chrBin,c(chr,binID),sep = "_")

# Preprocessing ------
## Matrix is symmetric. 
sparseM <- readMM('allGenesMatrix.mtx')

geneList <- unique(readChromDF$geneName)
interestingGenes <- unique(readChromDF %>% 
                             filter(geneType %in% c("protein_coding","lncRNA")) %>% 
                             .$geneName)

## Need  to account for distance dependent bias
## Next two lines because I think i made sparse mat with ALL genes but read chrom DF with just prot coding + lncRNA
genesToChoose <- which(geneList %in% interestingGenes)[1:1000]
allGeneNames <- geneList[genesToChoose]

newM <- sparseM[genesToChoose,genesToChoose]
reduced_binnedGeneDF <- inner_join(data.frame(geneName = allGeneNames),binnedGeneDF)
group = reduced_binnedGeneDF$gID

for(i in 1:max(reduced_binnedGeneDF$gID)){
  if(length(which(group == 1)) > 1){
    group = sort(c(i, group[group != i]))
    ix = which(reduced_binnedGeneDF$gID == i)
    B = Matrix(0,length(group),length(group))
    B[i,i] = sum(newM[ix,ix])
    B[i,-i] = colSums(newM[ix,-ix])
    B[-i,i] = rowSums(newM[-ix,ix])
    B[-i,-i] = newM[-ix,-ix]
    newM = B
    newDF = rbind(reduced_binnedGeneDF[reduced_binnedGeneDF$gID==i,][1,],
                  reduced_binnedGeneDF[reduced_binnedGeneDF$gID!=i,])
    df <- newDF[order(newDF[,2]),]
  }
}



#### Try
A = Matrix(rbinom(100,2,0.1),nrow = 10)
df <- data.frame("OG" = LETTERS[1:10], "New" = sort(sample(3,10,T)))
group = df$New

for(i in 1:max(df$New)){
  group = sort(c(i, group[group != i]))
  ix = which(df$New == i)
  B = Matrix(0,length(group),length(group))
  B[i,i] = sum(A[ix,ix])
  B[i,-i] = colSums(A[ix,-ix])
  B[-i,i] = rowSums(A[-ix,ix])
  B[-i,-i] = A[-ix,-ix]
  A = B
  newDF = rbind(df[df$New==i,][1,],df[df$New!=i,])
  df <- newDF[order(newDF[,2]),]
}



