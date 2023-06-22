#!/bin/R
# By Anoushka Joglekar 04.2023
# V2
## Reads in a sparse matrix containing pairwise gene x gene interactions
## from reads spanning multiple genes. 
## First, it constructs a network diagram, then it identifies gene x gene
## hubs. At no point am I considering the distance to the closest gene.
## Once appropriate hubs are identified, it plots the SE for those hub genes
## At present this is for GM12878 data only, but can be expanded to other cell lines / tissue

## Now that I'm thinking about it a tensor would be more appropriate. Otherwise
## we are not really taking advantage of the Pore-C nature other than knowing the cardinality per read.
## but let's start with baby steps.

## Also making sure to get rid of spurious connections

# Library loading ------
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
library(irlba)
library(leiden)
library(ComplexHeatmap)
library(uwot)

# Reading input ------
scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

SE_df <- read.table('../2023_03_06_GM12878_cellularFractionData/v0.evaluateSplicingEfficiency/SE_AvgOverReps_NuclearSE_withStatus', 
                    header = T) %>% dplyr::select(gene_ID,Nuclear,nucStatus)
readChromDF <- fread('../v0_poreC_explore/readChromunities_wClosestGene.gz')
colnames(readChromDF) <- c("ReadID","Cardinality","gene_ID","geneType","geneName","distToGene")

# Preprocessing ------
## Matrix is symmetric. Let's summarize the data
sparseM <- readMM('allGenesMatrix.mtx')
## Set diagonal elements to zero
diag(sparseM) <- 0

rS <- rowSums(sparseM)
threshold <- quantile(rS,0.25)
geneList <- unique(readChromDF$geneName)
protCoding <- unique(readChromDF %>% filter(geneType == "protein_coding") %>% .$geneName)


genesToChoose <- which(rS > threshold & geneList %in% protCoding)
newGeneNames <- geneList[genesToChoose]

newM <- sparseM[genesToChoose,genesToChoose]
colnames(newM) <- rownames(newM) <- newGeneNames

g <- graph_from_adjacency_matrix(newM, mode = "undirected")
g <- simplify(g,remove.multiple = F, remove.loops = T)

## getting imp hubs -------
deg <- degree(g, mode="all")
V(g)$size <- deg*0.001

hist(deg,breaks = 250,xlim = c(0,2000))
boxplot(deg,ylim = c(0,1000),notch = T, outline = F)

hub_threshold <- mean(deg) + 2 * sd(deg)
hub_threshold <- 10
hubs <- which(deg >= hub_threshold)

hubGeneNames <- newGeneNames[hubs]
hubGraph <- induced.subgraph(g,hubs)
hubNeighbors <- neighborhood(g,nodes = hubs)
names(hubNeighbors) <- names(hubs)

## can we identify communities or cliques? ----------
communities <- cluster_louvain(hubGraph)
membership <- communities$membership

## cliques
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
degHG <- degree(hubGraph)
#V(hubGraph)$size <- 8 #degHG*0.005

V(hubGraph)$size <- ifelse(degHG <= 80, 5, degHG * 0.01)
hubGraph <- delete_edges(hubGraph, E(hubGraph))
l = layout_with_fr(hubGraph)
plot(hubGraph, edge.arrow.size=.4, vertex.label = group_ids, vertex.label.cex = 0.5,
     vertex.color = colors[membership], vertex.label.color = colors[membership], 
     vertex.frame.color=NA,vertex.label.dist = 0.4)

## One big cluster and several very small ones. Let's plot the others individually ----
# First get the vertices in each component
mG <- NULL
for (i in unique(membership)) {
  vertices <- hubGeneNames[membership == i]
  mG[[i]] <- data.frame(gene_Name = vertices, cluster = i)
}
memberGenes <- do.call('rbind',mG)

maxHub <- hub_sizes %>% slice_max(size,n = 1) %>% .$hub_label
extraneousClusters <- memberGenes %>% filter(cluster != maxHub)
#deg[match(extraneousClusters$gene_Name,newGeneNames)]
n = neighborhood(g, nodes = extraneousClusters$gene_Name)

V(g)$size <- 0.001

layout(matrix(1:9, nrow = 3, ncol = 3))
for(i in 1:length(n)){
  s <- induced.subgraph(g,unique(as.vector(unlist(n[[i]]))))
  V(s)$color <- ifelse(V(s)$name == extraneousClusters$gene_Name[i], "red","black")
  #plot(s,label.cex = 0.1,vertex.size = 0.01,vertex.label.color = "black",vertex.label.dist = 0.7,vertex.frame.color=NA)
  plot(s,vertex.label = NA,vertex.size = 5,vertex.frame.color=NA)
}
dev.off()

## all together:
s <- induced.subgraph(g,unique(as.vector(unlist(n))))
V(s)$color <- ifelse(V(s)$name %in% extraneousClusters$gene_Name, "red","black")
l <- layout_with_fr(s)
plot(s,vertex.label = NA,vertex.size = 5,vertex.frame.color=NA)


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
hubMemberGenes <- do.call('rbind',hmG) ## have not looked at overlaps

## tying back to SE ---------
big_hubs <- hub_sizes %>% filter(size >=5)
annot <- read.table('../v0_poreC_explore/geneCoordsAndNames_sorted.bed',sep = "\t")[,c(5,7)]
annot$V5 <- gsub(";$", "", annot$V5)
annot$V7 <- gsub(";$", "", annot$V7)
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
  geom_point(data = SE_df_wClustIDs %>% filter(hubStatus == "Hub"), color = "red", alpha = 0.2) +
  ylab("SE per gene by cluster")

## tying back to SE part 2 i.e. by hubs ---------
SE_df_wHubNames <- inner_join(hubMemberGenes,SE_df)
SE_df_wHubNames <- SE_df_wHubNames%>% group_by(hubGeneName) %>%
  summarize(numGenes = n(),meanSE = median(Nuclear)) %>%
  ungroup() %>% mutate(binned_numGenes = cut(numGenes,breaks = seq(0,1100,by = 100) )) %>%
  group_by(binned_numGenes) %>%
  add_count(name = "numHubs")


#SE_df_wHubNames$binned_numGenes <- as.factor(SE_df_wHubNames$binned_numGenes)
ggplot(SE_df_wHubNames, #%>% filter(numHubs >= 5),
       aes(x = binned_numGenes,y = meanSE)) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75)) +
  ylim(0.7,1) + ggtitle("Distribution of splicing efficiency across hubs") +
  geom_signif(comparisons = list(c("(0,100]","(100,200]")), y_position = c(0.97))



## trying dim reduction algorithms like PCA first
dS <- diag(sparseM)
diag(sparseM) <- 0

rS <- rowSums(sparseM)
#threshold <- quantile(rS,0.25)
geneList <- unique(readChromDF$geneName)
protCoding <- unique(readChromDF %>% filter(geneType == "protein_coding") %>% .$geneName)


genesToChoose <- which(geneList %in% protCoding)
newGeneNames <- geneList[genesToChoose]

newM <- sparseM[genesToChoose,genesToChoose]
colnames(newM) <- rownames(newM) <- newGeneNames
nonZero <- union(which(colSums(newM) != 0),which(rowSums(newM) != 0))
newGeneSet <- newGeneNames[nonZero]
newNewM <- newM[newGeneSet,newGeneSet]
norm_newM <- sweep(newNewM,2,colSums(newNewM),`/`)

set.seed(1)
randomGenes <- intersect(chr18$geneName,newGeneNames) #sample(newGeneNames,300,replace = F)

subNewM <- newM[randomGenes,randomGenes]
nonZero <- union(which(colSums(subNewM) != 0),which(rowSums(subNewM) != 0))
newGeneSet <- randomGenes[nonZero]
subsub <- subNewM[newGeneSet,newGeneSet]
norm_newM <- sweep(subsub,2,colSums(subsub),`/`)

norm_dense <- as.matrix(norm_newM)

colFunc <- circlize::colorRamp2(c(-0.25, 0, 0.25), c("blue","white", "red"))
corNorm_dense <- cor(norm_dense)
Heatmap(corNorm_dense,
        cluster_rows = F, cluster_columns = F,
        show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F,
        name = "chr18",col = colFunc)

pca_corNorm <- prcomp_irlba(corNorm_dense, n=20)
sdev <- pca_corNorm$sdev
plot(1:length(sdev), sdev^2, type = "b", 
     xlab = "Principal Component", ylab = "Variance Explained", pch = 16)

scores <- pca_corNorm$x
plot(scores[,1],scores[,2],pch = 16)
k <- 4

kmeans_result <- kmeans(scores, centers = k)
cluster_labels <- kmeans_result$cluster
plot(scores[, 1], scores[, 2], col = cluster_labels, pch = 16, main = "PCA - K-means Clustering")   

# chrIx <- left_join(data.frame(geneName = rownames(corNorm_dense)),geneCoords) %>% 
#   dplyr::select(chr,geneName) %>% group_by(chr) %>% mutate(group_id = cur_group_id())
# 
# plot(scores[, 1], scores[, 2], col = chrIx$group_id, pch = 16, main = "PCA - Chr ID")

## umap on data
umap_corNorm <- umap(scores, 8,nn_method = "annoy")
k <- 5
kmeans_result <- kmeans(umap_corNorm, centers = k)
cluster_labels <- data.frame(geneName=rownames(corNorm_dense), clusterID=kmeans_result$cluster)
plot(umap_corNorm[, 1], umap_corNorm[, 2], pch = 16, main = "UMAP", col = cluster_labels$clusterID)   
#plot(umap_corNorm[, 1], umap_corNorm[, 2], pch = 16, main = "UMAP", col = chrIx$group_id)   

# ## distance matrix from pca
# distance_matrix <- dist(scores)
# leiden_result <- leiden(distance_matrix)
# louvain_result <- cluster_louvain(distance_matrix)


## get graph
splitDF <- split.data.frame(cluster_labels,cluster_labels$clusterID)
clust2_genes <- splitDF$`2`$geneName

filteredCorMat <- newM[clust2_genes,clust2_genes]
#filteredCorMat[filteredCorMat <= 0] <- 0
g <- graph_from_adjacency_matrix(filteredCorMat, mode = "undirected",diag = F)
g <- simplify(g,remove.multiple = F, remove.loops = T)

deg <- degree(g, mode="all")
V(g)$size <- deg*0.01

hist(deg,breaks = 200,xlim = c(0,100))

hub_threshold <- 10
hubs <- which(deg >= hub_threshold)

hubGeneNames <- clust8_genes[hubs]
hubGraph <- induced.subgraph(g,hubs)
hubNeighbors <- neighborhood(g,nodes = hubs)
names(hubNeighbors) <- names(hubs)

## can we identify communities or cliques? ----------
communities <- cluster_louvain(hubGraph)
membership <- data.frame(geneID=names(hubs),ogCluster="clust1",memberID = communities$membership)



colors <- viridis_pal(option = "viridis")(length(unique(membership$memberID)))
degHG <- degree(hubGraph)
V(hubGraph)$size <- 8 #degHG*0.005
#hubGraph <- delete_edges(hubGraph, E(hubGraph))
l = layout_with_fr(hubGraph)
plot(hubGraph, edge.arrow.size=.4, vertex.label.cex = 0.8,
     #layout = l, #edge.color = white,vertex.label = membership,
     vertex.color = colors[membership$memberID], vertex.label.color = colors[membership$memberID], 
     vertex.frame.color=NA,vertex.label.dist = 0.4)


## Combine the graph with the complex heatmap to get annotations
cluster_labels$clusterID <- as.factor(cluster_labels$clusterID)
hubAnnot <- left_join(data.frame(geneID=rownames(corNorm_dense)),membership) %>% 
  select(-ogCluster) %>% replace(is.na(.),"None")
hM_annotation <- HeatmapAnnotation(umapClusters = cluster_labels$clusterID,
                                   hubAnnot = hubAnnot$memberID)
Heatmap(corNorm_dense,
        cluster_rows = F, cluster_columns = F,
        show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F,
        top_annotation = hM_annotation,
        name = "chr18", col =colFunc)

## for HHW
clust17 <- memberGenes[memberGenes$cluster == 17,]
library(enrichR)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched <- enrichr(clust17$gene_Name,dbs)
enrichR::plotEnrich(enriched$GO_Biological_Process_2015)
## get neighborhood of all genes in hub17

clust17_sparse <- newM[clust17$gene_Name,clust17$gene_Name]

s <- induced.subgraph(g,clust17$gene_Name)
clusterlouvain <- cluster_louvain(s)

group_ids <- clusterlouvain$membership
duplicated_indices <- duplicated(group_ids)
group_ids[duplicated_indices] <- NA
rm(duplicated_indices)

l <- layout_nicely(s)
#s <- delete_edges(s, E(s))
plot(s,vertex.size = 8,vertex.frame.color=NA,vertex.label.cex = 2,vertex.label = NA,
     #vertex.label = group_ids,
     layout = l,
     #vertex.label.color=viridis(option = "D",n = length(unique(group_ids)))[clusterlouvain$membership],
     vertex.color=colorRampPalette(brewer.pal(5,"BuGn"))(length(unique(group_ids)))[clusterlouvain$membership])

> which(clust17_sparse == max(clust17_sparse),arr.ind = T)
# row col
# RBM45  111 110
# PDE11A 110 111

n <- neighborhood(g,nodes = "RBM45")
s <- induced.subgraph(g,"RBM45")
V(s)$color <- ifelse(V(s)$name == "RBM45", "red","black")
#plot(s,label.cex = 0.1,vertex.size = 0.01,vertex.label.color = "black",vertex.label.dist = 0.7,vertex.frame.color=NA)
plot(s,vertex.label = NA,vertex.size = 5,vertex.frame.color=NA)

head(memberGenes %>% group_by(cluster) %>% n_distinct())


## power law function:
geneCoords <- read.table('../v0_poreC_explore/geneCoordsAndNames_sorted.bed')
colnames(geneCoords) <- c("chr","start","end","strand","geneID",
                          "s1","type","s2","geneName","s3")
geneCoords <- geneCoords %>% dplyr::select(-c(s1,s2,s3))

chromSizes <- read.table('../v0_poreC_explore/hg38.chromSizes')
colnames(chromSizes) <- c("chr","size")

chr18 <- geneCoords %>% filter(chr %in% c("chr18"))

chr18 <- geneCoords %>% filter(chr == "chr18") %>% 
  mutate_at(c("start","end"),as.integer) %>%
  mutate( coordStart = cut( start, breaks = seq(1,chromSizes[18,2],10000) ),
          coordEnd = cut( end, breaks = seq(1,chromSizes[18,2],10000) ))


chr17 <- geneCoords %>% filter(chr == "chr17") %>% 
  mutate_at(c("start","end"),as.integer) %>%
  mutate( coordStart = cut( start, breaks = seq(1,chromSizes[17,2],10000) ),
          coordEnd = cut( end, breaks = seq(1,chromSizes[17,2],10000) ))


pInt <- function(x,y){
  p = dist(x,y)
}






