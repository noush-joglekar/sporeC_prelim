# !/bin/R
# R version 4.1.2
# Automated differential analysis result interpretation + visualization w/ gTrack

## Load packages
library(dplyr)
library(tidyr)
library(igraph)

library(gTrack)
library(rtracklayer)

library(ComplexHeatmap)
library(circlize)

args <- commandArgs(trailing = TRUE)

# Setup -----
# args <- c('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/',
#           'differentialAnalysis/gm_vs_hcc/','GM12878','HCC1954','3','chr11','differentialAnalysis/Plots_GM_vs_HCC/')

workingDir <- args[1]
setwd(workingDir)
daDir <- args[2]

ct1 <- args[3]
ct2 <- args[4]

card <- args[5]
chrom <- args[6]

plotDir <- args[7]
if(!dir.exists(plotDir)){dir.create(plotDir)}

## Cell type 1 ----------
impDiffReads <- read.table(file.path(workingDir,daDir,
                                     paste0('/diffContactReads_card',card,'_',chrom,'_',ct1,'.tab.gz')), header = T)
commAnnot <- read.table(paste0(daDir,'binCommunities_',ct1,'_card',card,'_',chrom,'.csv'),
                        sep=",", header = T)  %>% 
  separate_rows(Nodes,sep = ",") %>% dplyr::rename(Cluster = Community.ID, binID = Nodes)
impReads_wCluster <- right_join(impDiffReads,commAnnot)

impReads_reformatted <- as.data.frame(impReads_wCluster) %>% dplyr::select(chr,binStart,binEnd,binID,Cluster,geneName,readID)
colnames(impReads_reformatted) <- c("seqnames","start","end","binID","Cluster","Gene","readID")

ct1_gr <- dt2gr(impReads_reformatted)

ic_gr <- list()
for (i in sort(unique(commAnnot$Cluster))){
  ic_gr[i+1] <- ct1_gr %Q% (Cluster == i)
}

my_string = c()
for(i in sort(unique(commAnnot$Cluster))){
  my_string = c(my_string,paste0("gTrack(split(ic_gr[",i+1,"][[1]], ic_gr[",i+1,"][[1]]$readID) %>% unname, height = 10, name = '",
                                 paste0('C',i),"')"))
}

my_string_ct1 = gsub("\"","",paste(my_string,collapse = ","))
toEval_ct1 <- paste0("plot(c(",my_string_ct1,"))")

pdf(paste0(plotDir,ct1,'_clusters_vs',ct2,'_card',card,'_',chrom,'.pdf'),12,10,useDingbats = F)
eval(parse(text=toEval_ct1))
dev.off()


## Cell type 2 ------
impDiffReads_ct2 <- read.table(file.path(workingDir,daDir,
                                         paste0('/diffContactReads_card',card,'_',chrom,'_',ct2,'.tab.gz')), header = T)
commAnnot_ct2 <- read.table(paste0(daDir,'binCommunities_',ct2,'_card',card,'_',chrom,'.csv'),
                           sep=",", header = T)  %>% 
  separate_rows(Nodes,sep = ",") %>% dplyr::rename(Cluster = Community.ID, binID = Nodes)
impReads_wCluster_ct2 <- right_join(impDiffReads_ct2,commAnnot_ct2)

impReads_reformatted_ct2 <- as.data.frame(impReads_wCluster_ct2) %>% dplyr::select(chr,binStart,binEnd,binID,Cluster,geneName,readID)
colnames(impReads_reformatted_ct2) <- c("seqnames","start","end","binID","Cluster","Gene","readID")

ct2_gr <- dt2gr(impReads_reformatted_ct2)

ic_gr_ct2 <- list()
for (i in sort(unique(commAnnot_ct2$Cluster))){
  ic_gr_ct2[i+1] <- ct2_gr %Q% (Cluster == i)
}

my_string = c()
for(i in sort(unique(commAnnot_ct2$Cluster))){
  my_string = c(my_string,paste0("gTrack(split(ic_gr_ct2[",i+1,"][[1]], ic_gr_ct2[",i+1,"][[1]]$readID) %>% unname, height = 10, name = '",
                                 paste0('C',i),"')"))
}

my_string_ct2 = gsub("\"","",paste(my_string,collapse = ","))
toEval_ct2 <- paste0("plot(c(",my_string_ct2,"))")

pdf(paste0(plotDir,ct2,'_clusters_vs',ct1,'_card',card,'_',chrom,'.pdf'),12,10,useDingbats = F)
eval(parse(text=toEval_ct2))
dev.off()


## Plot both full GRanges together -----
pdf(paste0(plotDir,'both_',ct1,'_',ct2,'_card',card,'_',chrom,'_fullGR.pdf'),15,12,useDingbats = F)
plot(c(gTrack(split(ct1_gr, ct1_gr$readID) %>% unname, height = 10, name = ct1),
       gTrack(split(ct2_gr, ct2_gr$readID) %>% unname, height = 10, name = ct2)))
dev.off()


## Get corresponding clusters (not sure why) ------

jacMat <- matrix(0, length(ic_gr), length(ic_gr_ct2))
rownames(jacMat) <- paste0("C",(1:length(ic_gr)-1))
colnames(jacMat) <- paste0("C",(1:length(ic_gr_ct2)-1))

for (c1 in 1:length(ic_gr)){
  gl1 = ic_gr[[c1]]$Gene
  for (c2 in 1:length(ic_gr_ct2)){
    gl2 = ic_gr_ct2[[c2]]$Gene
    li = length(intersect(gl1,gl2))
    lu = length(union(gl1,gl2))
    jacMat[c1,c2] <- (li / lu)
  }
}

col_fun <- circlize::colorRamp2(c(0,0.5,1),colors = c("white","pink","maroon"))
Heatmap(jacMat,name = chrom, col = col_fun,row_title = ct1,column_title = ct2, border = "grey50",
        show_column_dend = F, show_row_dend = F)




