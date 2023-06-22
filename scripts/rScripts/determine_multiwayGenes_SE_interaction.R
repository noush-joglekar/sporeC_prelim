#!/bin/R
# By Anoushka Joglekar 04.2023
## Reads in a sparse matrix containing pairwise gene x gene interactions
## from reads spanning multiple genes. 
## First, it constructs a network diagram, then it identifies gene x gene
## hubs. At no point am I considering the distance to the closest gene.
## Once appropriate hubs are identified, it plots the SE for those hub genes
## At present this is for GM12878 data only, but can be expanded to other cell lines / tissue

## Now that I'm thinking about it a tensor would be more appropriate. Otherwise
## we are not really taking advantage of the Pore-C nature other than knowing the cardinality per read.
## but let's start with baby steps.

library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(parallel)
library(igraph)
library(ggplot2)
library(ggsignif)
library(viridis)
library(gridExtra)
library(MASS)

scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

SE_df <- read.table('../2023_03_06_GM12878_cellularFractionData/v0.evaluateSplicingEfficiency/SE_AvgOverReps_NuclearSE_withStatus', 
                    header = T) %>% select(gene_ID,Nuclear,nucStatus)
readChromDF <- fread('../v0_poreC_explore/readChromunities_wClosestGene.gz')
colnames(readChromDF) <- c("ReadID","Cardinality","gene_ID","geneName","distToGene")

## Network diagram attempt -------
geneList <- unique(readChromDF$geneName)
sparseM <- readMM('allGenesMatrix.mtx')
colnames(sparseM) <- rownames(sparseM) <- geneList
nonZero <- union(which(colSums(sparseM) != 0),which(rowSums(sparseM) != 0))
newGeneNames <- geneList[nonZero]

newM <- sparseM[nonZero,nonZero]
colnames(newM) <- rownames(newM) <- newGeneNames
rm(sparseM)
g <- graph_from_adjacency_matrix(newM, mode = "undirected")
g <- simplify(g,remove.multiple = F, remove.loops = T)

deg <- degree(g, mode="all")
V(g)$size <- deg*0.01

## I think the next part may be a bit much for plotting: -----
# Create force-directed layout
l <- layout_with_fr(g)
# Spherical layout
#l <- layout_on_sphere(g)
E(g)$edge.color <- NA
plot(g, edge.arrow.size=.4,vertex.label=NA,
     layout = l,
     vertex.color = "blue", vertex.frame.color=NA)

## getting imp hubs -------
hist(deg,breaks = 500,xlim = c(0,2000))
boxplot(deg,ylim = c(0,1500),notch = T, outline = F)
# fitDeg <- fitdistr(deg, "negative binomial")
# > fitDeg$estimate[1]
# size 
# 0.6696002 
# > fitDeg$estimate[2]
# mu 
# 432.5528 

#hub_threshold <- mean(deg) + 2 * sd(deg)
hub_threshold <- quantile(deg, 0.95)
# hub_min <- quantile(deg, 0.5)
# hub_max <- quantile(deg, 0.75)
# hubs <- which(deg > hub_min & deg < hub_max)
hubs <- which(deg >= hub_threshold)

hubGeneNames <- newGeneNames[hubs]
hubGraph <- induced.subgraph(g,hubs)
hubNeighbors <- neighborhood(g,nodes = hubs)
names(hubNeighbors) <- names(hubs)

## all together (might be too much again) ------
hubsWithNeighbors <- induced.subgraph(g,unique(as.vector(unlist(hubNeighbors))))
V(hubsWithNeighbors)$color <- ifelse(V(hubsWithNeighbors)$name %in% hubGeneNames, "red","black")
l <- layout_with_fr(hubsWithNeighbors)
hubsWithNeighbors <- delete_edges(hubsWithNeighbors, E(hubsWithNeighbors))
plot(hubsWithNeighbors,vertex.label = NA,vertex.frame.color=NA,layout = l)

# compute the connected components of the graph -----
comps <- components(hubGraph)

# get the membership vector
membership <- comps$membership

hub_sizes <- as.data.frame(table(membership))
colnames(hub_sizes) <- c("hub_label", "size")

# get the vertices in each component
mG <- NULL
for (i in unique(membership)) {
  vertices <- hubGeneNames[membership == i]
  mG[[i]] <- data.frame(gene_Name = vertices, cluster = i)
}
memberGenes <- do.call('rbind',mG)

## Label and color plot by membership group ----
group_ids <- membership
duplicated_indices <- duplicated(group_ids)
group_ids[duplicated_indices] <- NA
rm(duplicated_indices)

colors <- viridis_pal(option = "viridis")(length(unique(membership)))

plot(hubGraph, edge.arrow.size=.4, vertex.label = group_ids, vertex.label.cex = 1,
     vertex.color = colors[membership], vertex.label.color = colors[membership], 
     vertex.frame.color=NA,vertex.label.dist = 0.7)

## One big cluster and several very small ones. Let's plot the others individually. 
maxHub <- hub_sizes %>% slice_max(size,n = 1) %>% .$hub_label
extraneousClusters <- memberGenes %>% filter(cluster != maxHub)
#deg[match(extraneousClusters$gene_Name,newGeneNames)]
n = neighborhood(g, nodes = extraneousClusters$gene_Name)

layout(matrix(1:6, nrow = 2, ncol = 3))
for(i in 1:length(n)){
  s <- induced.subgraph(g,unique(as.vector(unlist(n[[i]]))))
  V(s)$color <- ifelse(V(s)$name == extraneousClusters$gene_Name[i], "red","black")
  plot(s,label.cex = 0.6)
}
dev.off()
## all together:
s <- induced.subgraph(g,unique(as.vector(unlist(n))))
V(s)$color <- ifelse(V(s)$name %in% extraneousClusters$gene_Name, "red","black")
l <- layout_with_fr(s)
plot(s,vertex.label = NA,layout = l)

## can we identify communities or cliques and re-plot? ----------
## works a bit better
igraph::min_cut(hubGraph, value.only = FALSE)
communities <- cluster_louvain(hubGraph)
membership <- communities$membership

## OR
cliques <- max_cliques(hubGraph, min = 10, max = 30)
CEB <- cluster_edge_betweenness(hubGraph)
membership <- CEB$membership

hub_sizes <- as.data.frame(table(membership))
colnames(hub_sizes) <- c("hub_label", "size")

group_ids <- membership
duplicated_indices <- duplicated(group_ids)
group_ids[duplicated_indices] <- NA
rm(duplicated_indices)

colors <- viridis_pal(option = "viridis")(length(unique(membership)))

V(hubGraph)$size <- 0.01
hubGraph <- delete_edges(hubGraph, E(hubGraph))
plot(hubGraph, edge.arrow.size=.4, vertex.label = group_ids, vertex.label.cex = 0.5,
     vertex.color = colors[membership], vertex.label.color = colors[membership], 
     vertex.frame.color=NA,vertex.label.dist = 0.4)

# get the vertices in each clique ------
mG <- NULL
for (i in unique(membership)) {
  vertices <- hubGeneNames[membership == i]
  connectedGenes <- unique(unlist(unname(sapply(vertices, function(v) names(hubNeighbors[[v]])))))
  mG[[i]] <- data.frame(gene_Name = connectedGenes, cluster = i) %>% 
    mutate(hubStatus = case_when(gene_Name %in% vertices ~ "Hub", TRUE ~ "NonHub"))
}
memberGenes <- do.call('rbind',mG)

## get all the genes in a hub -----
hmG <- NULL
for (hg in hubGeneNames) {
  connectedGenes <- names(hubNeighbors[[hg]])
  hmG[[hg]] <- data.frame(gene_Name = connectedGenes, hubGeneName = hg) %>% 
    mutate(hubStatus = case_when(gene_Name == hg ~ "Hub", TRUE ~ "NonHub"))
}
hubMemberGenes <- do.call('rbind',hmG)

## tying back to SE ---------
big_hubs <- hub_sizes %>% filter(size >=10)
annot <- read.table('../v0_poreC_explore/geneCoordsAndNames_sorted.bed',sep = "\t")[,c(5,6)]
annot$V5 <- gsub(";$", "", annot$V5)
annot$V6 <- gsub(";$", "", annot$V6)
colnames(annot) <- c("gene_ID","gene_Name")

SE_df <- left_join(SE_df,annot)

SE_df_wClustIDs <- inner_join(memberGenes %>% filter(cluster %in% big_hubs$hub_label),SE_df)
SE_df_wClustIDs$cluster <- as.factor(SE_df_wClustIDs$cluster)

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
}

ggplot(SE_df_wClustIDs,aes(x = cluster,y = Nuclear)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75)) +
  geom_point(data = SE_df_wClustIDs %>% filter(hubStatus == "Hub"), color = "red", alpha = 0.2)

## tying back to SE part 2 i.e. by hubs ---------
SE_df_wHubNames <- inner_join(hubMemberGenes,SE_df)
SE_df_wHubNames <- SE_df_wHubNames%>% group_by(hubGeneName) %>%
  summarize(numGenes = n(),meanSE = median(Nuclear)) %>%
  ungroup() %>% group_by(numGenes) %>%
  add_count(name = "numHubs")


give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
}

SE_df_wHubNames$numGenes <- as.factor(SE_df_wHubNames$numGenes)
ggplot(SE_df_wHubNames %>% filter(numHubs >= 5),
       aes(x = numGenes,y = meanSE)) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75)) +
  ylim(0.7,1) + ggtitle("Distribution of splicing efficiency across hubs") +
  geom_signif(comparisons = list(c("16","26"),c("21","26")), y_position = c(0.97,0.95))

wilcox.test(SE_df_wHubNames %>% filter(numGenes == 6) %>% .$meanSE,
        SE_df_wHubNames %>% filter(numGenes == 7) %>% .$meanSE)


