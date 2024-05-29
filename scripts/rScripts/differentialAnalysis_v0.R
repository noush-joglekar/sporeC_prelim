#!/bin/R
# R version 4.1.2
# Differential analysis result interpretation

## Load packages
library(dplyr)
library(tidyr)
library(igraph)

library(gTrack)

# Setup -----
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/'
setwd(workingDir)

ct1 <- 'HCC1954'
ct2 <- 'HG002'

impDiffReads <- read.table('differentialAnalysis/gm_vs_hcc//diffContactReads_card3_chr11_HCC1954.tab.gz', header = T)
commAnnot <- read.table('differentialAnalysis/gm_vs_hcc//binCommunities_HCC1954_card3_chr11.csv',sep=",", header = T)  %>% 
  separate_rows(Nodes,sep = ",") %>% rename(Cluster = Community.ID, binID = Nodes)
impReads_wCluster <- right_join(impDiffReads,commAnnot)

graph <- read_graph('differentialAnalysis/hcc_vs_hg/weightedLineGraph_HCC1954_readByBin_card3_chr8.graphml', format = "graphml")


impReads_reformatted <- as.data.frame(impReads_wCluster) %>% dplyr::select(chr,binStart,binEnd,binID,Cluster,geneName,readID)
colnames(impReads_reformatted) <- c("seqnames","start","end","binID","Cluster","Gene","readID")

hcc_gr <- dt2gr(impReads_reformatted)

ic_gr <- list()
for (i in sort(unique(commAnnot$Cluster))){
  ic_gr[i+1] <- hcc_gr %Q% (Cluster == i)
}

view_range = GRanges('chr11:500001-135000000')


my_string = c()
for(i in sort(unique(commAnnot$Cluster))){
  my_string = c(my_string,paste0("gTrack(split(ic_gr[",i+1,"][[1]], ic_gr[",i+1,"][[1]]$readID) %>% unname, height = 10, name = '",paste0('C',i),"')"))
}

my_string = gsub("\"","",paste(my_string,collapse = ","))
my_string
plot(c(my_string),view_range)
plot(c(#gTrack(split(ic_gr[0][[1]], ic_gr[0][[1]]$readID) %>% unname, height = 10, name = 'C0'),
     gTrack(split(ic_gr[1][[1]], ic_gr[1][[1]]$readID) %>% unname, height = 10, name = 'C1'),
     gTrack(split(ic_gr[2][[1]], ic_gr[2][[1]]$readID) %>% unname, height = 10, name = 'C2'),
     gTrack(split(ic_gr[3][[1]], ic_gr[3][[1]]$readID) %>% unname, height = 10, name = 'C3'),
     gTrack(split(ic_gr[4][[1]], ic_gr[4][[1]]$readID) %>% unname, height = 10, name = 'C4')))

#### -----
# G = plot(c(gTrack(split(ic_gr0, ic_gr0$readID) %>% unname, height = 10, name = 'C0'),
#   gTrack(split(ic_gr1, ic_gr1$readID) %>% unname, height = 10, name = 'C1'),
#   gTrack(split(ic_gr2, ic_gr2$readID) %>% unname, height = 10, name = 'C2'),
#   gTrack(split(ic_gr3, ic_gr3$readID) %>% unname, height = 10, name = 'C3'),
#   gTrack(split(ic_gr4, ic_gr4$readID) %>% unname, height = 10, name = 'C4'),
#   gTrack(split(ic_gr5, ic_gr5$readID) %>% unname, height = 10, name = 'C5'),
#   gTrack(split(ic_gr6, ic_gr6$readID) %>% unname, height = 10, name = 'C6'),
#   gTrack(split(ic_gr7, ic_gr7$readID) %>% unname, height = 10, name = 'C7'),
#   gTrack(split(ic_gr8, ic_gr8$readID) %>% unname, height = 10, name = 'C8'),
#   gTrack(split(ic_gr9, ic_gr9$readID) %>% unname, height = 10, name = 'C9'),
#   gTrack(split(ic_gr10, ic_gr10$readID) %>% unname, height = 10, name = 'C10'),
#   gTrack(split(ic_gr11, ic_gr11$readID) %>% unname, height = 10, name = 'C11'),
#   gTrack(split(ic_gr12, ic_gr12$readID) %>% unname, height = 10, name = 'C12'),
#   gTrack(split(ic_gr13, ic_gr13$readID) %>% unname, height = 10, name = 'C13')
#   ),view_range)


### HG002 -----
impDiffReads_hg <- read.table('differentialAnalysis/hcc_vs_hg/diffContactReads_card3_chr11_HG002.tab.gz', header = T)
commAnnot_hg <- read.table('differentialAnalysis/hcc_vs_hg/binCommunities_HG002_card3_chr11.csv',sep=",", header = T)  %>% 
  separate_rows(Nodes,sep = ",") %>% rename(Cluster = Community.ID, binID = Nodes)
impReads_wCluster_hg <- right_join(impDiffReads_hg,commAnnot_hg)

impReads_reformatted_hg <- as.data.frame(impReads_wCluster_hg) %>% dplyr::select(chr,binStart,binEnd,binID,Cluster,geneName,readID)
colnames(impReads_reformatted_hg) <- c("seqnames","start","end","binID","Cluster","Gene","readID")

hg_gr <- dt2gr(impReads_reformatted_hg)

ic_gr_hg <- list()
for (i in sort(unique(commAnnot_hg$Cluster))){
  ic_gr_hg[i+1] <- hg_gr %Q% (Cluster == i)
}

my_string = c()
for(i in sort(unique(commAnnot_hg$Cluster))){
  my_string = c(my_string,paste0("gTrack(split(ic_gr_hg[",i+1,"][[1]], ic_gr_hg[",i+1,"][[1]]$readID) %>% unname, height = 10, name = '",paste0('C',i),"')"))
}

my_string = gsub("\"","",paste(my_string,collapse = ","))
plot(c(my_string),view_range)

plot(c(#gTrack(split(ic_gr_hg[0][[1]], ic_gr_hg[0][[1]]$readID) %>% unname, height = 10, name = 'C0'),
       gTrack(split(ic_gr_hg[1][[1]], ic_gr_hg[1][[1]]$readID) %>% unname, height = 10, name = 'C1'),
       gTrack(split(ic_gr_hg[2][[1]], ic_gr_hg[2][[1]]$readID) %>% unname, height = 10, name = 'C2'),
       gTrack(split(ic_gr_hg[3][[1]], ic_gr_hg[3][[1]]$readID) %>% unname, height = 10, name = 'C3'),
       gTrack(split(ic_gr_hg[4][[1]], ic_gr_hg[4][[1]]$readID) %>% unname, height = 10, name = 'C4'),
       gTrack(split(ic_gr_hg[5][[1]], ic_gr_hg[5][[1]]$readID) %>% unname, height = 10, name = 'C5'),
       gTrack(split(ic_gr_hg[6][[1]], ic_gr_hg[6][[1]]$readID) %>% unname, height = 10, name = 'C6')))


### GM12878 -----
impDiffReads_gm <- read.table('differentialAnalysis/gm_vs_hcc/diffContactReads_card3_chr11_GM12878.tab.gz', header = T)
commAnnot_gm <- read.table('differentialAnalysis/gm_vs_hcc/binCommunities_GM12878_card3_chr11.csv',sep=",", header = T)  %>% 
  separate_rows(Nodes,sep = ",") %>% rename(Cluster = Community.ID, binID = Nodes)
impReads_wCluster_gm <- right_join(impDiffReads_gm,commAnnot_gm)

impReads_reformatted_gm <- as.data.frame(impReads_wCluster_gm) %>%
  dplyr::select(chr,binStart,binEnd,binID,Cluster,geneName,readID)
colnames(impReads_reformatted_gm) <- c("seqnames","start","end","binID","Cluster","Gene","readID")

gm_gr <- dt2gr(impReads_reformatted_gm)

ic_gr_gm <- list()
for (i in sort(unique(commAnnot_gm$Cluster))){
  ic_gr_gm[i+1] <- gm_gr %Q% (Cluster == i)
}

my_string = c()
for(i in sort(unique(commAnnot_gm$Cluster))){
  my_string = c(my_string,paste0("gTrack(split(ic_gr_gm[",i+1,"][[1]], ic_gr_gm[",i+1,"][[1]]$readID) %>% unname, height = 10, name = '",paste0('C',i),"')"))
}

my_string = gsub("\"","",paste(my_string,collapse = ","))
plot(c(my_string),view_range)

## BOTH:

plot(c(gTrack(split(hg_gr, hg_gr$readID) %>% unname, height = 10, name = 'HG002'),
       gTrack(split(gm_gr, gm_gr$readID) %>% unname, height = 10, name = 'GM12878'),
       gTrack(split(hcc_gr, hcc_gr$readID) %>% unname, height = 10, name = 'HCC1954')))

### SE for GM12878 ---------
# allReads_gm <- data.table::fread('NlaIII_GM12878_output_byChr/NlaIII_GM12878_chr11.gz', header = F)
nullAnnot_gm <- read.table('differentialAnalysis/gm_vs_hcc_test_unweightedbinCommunities_GM12878_card3_chr11.csv',
                           sep=",", header = T)  %>%
  separate_rows(Nodes,sep = ",") %>% rename(Cluster = Community.ID, binID = Nodes)

nullReads_wCluster_gm <- right_join(allReads_gm,nullAnnot_gm)

nullReads_reformatted_gm <- as.data.frame(nullReads_wCluster_gm) %>%
  dplyr::select(chr,binStart,binEnd,binID,Cluster,geneName,readID)
colnames(nullReads_reformatted_gm) <- c("seqnames","start","end","binID","Cluster","Gene","readID")

gm_gr_null <- dt2gr(nullReads_reformatted_gm)

ic_gr_gm_null <- list()
for (i in sort(unique(nullAnnot_gm$Cluster))){
  ic_gr_gm_null[i+1] <- gm_gr_null %Q% (Cluster == i)
}



SE_genes <- read.table('Gene_NucSE',sep = "\t")
colnames(SE_genes) <- c("GeneID","Gene","SE")

seList <- list()

for(i in sort(unique(commAnnot_gm$Cluster))){
  seList[[i+1]] = data.frame(SE = SE_genes %>% filter(Gene %in% ic_gr_gm[[i+1]]$Gene) %>% .$SE, Cluster = i)
}
seDF <- do.call('rbind',seList)
seDF$Type <- "weighted"

seNULL <- list()
for(i in sort(unique(nullAnnot_gm$Cluster))){
  seNULL[[i+1]] = data.frame(SE = SE_genes %>% filter(Gene %in% ic_gr_gm_null[[i+1]]$Gene) %>% .$SE, Cluster = i)
}
seNULL_DF <- do.call('rbind',seNULL)
seNULL_DF$Type <- "Null"

fullDF <- rbind(seDF,seNULL_DF)
fullDF$Cluster <- as.factor(fullDF$Cluster)

ggplot(fullDF, aes(x = Cluster, y = SE, color = Type)) +
  geom_boxplot() +ylim(0.85,1) +
  facet_wrap(~Type,scales = "free_x")

