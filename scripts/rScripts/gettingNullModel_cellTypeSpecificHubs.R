library(ComplexHeatmap)
library(circlize)

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
Heatmap(jacMat,name = "chr11", col = col_fun,row_title = "CT1",column_title = "CT2", border = "grey50",
        show_column_dend = F, show_row_dend = F)

ic_gr[[5]]
ic_gr_ct2[[9]]


## Splicing efficiency
SE_genes <- read.table('Gene_NucSE',sep = "\t")
colnames(SE_genes) <- c("GeneID","Gene","SE")

seList <- list()
for(i in sort(unique(commAnnot$Cluster))){
  seList[[i+1]] = data.frame(SE = SE_genes %>% filter(Gene %in% ic_gr[[i+1]]$Gene) %>% .$SE, Cluster = i)
}
seDF <- do.call('rbind',seList)
seDF$Type <- "Weighted"


## Null model:
gm <- fread('NlaIII_GM12878_output_byChr/NlaIII_GM12878_chr11.gz',nThread = 4)
gm <- gm[,c(1:4,10,11,12,13)]
colnames(gm) <- c("chr","start","end","seqnames","strand","geneID","bioType","Gene")
binnedGM <- gm %>% mutate(bin = cut_width(start, width = 500000, center = 250000)) %>% 
  group_by(bin) %>%
  mutate(binID = paste0("Bin",cur_group_id()))

set.seed(101)
SE_null <- list()
for (i in sort(unique(commAnnot$Cluster))){
  print(i)
  a <- as.integer(gsub("Bin","",unique(ic_gr[[i+1]]$binID)))
  b <- c()
  for(x in 1:length(a)){b <- c(b,a[x] - a[1])}
  minBin <- 1
  maxBin <- as.integer(gsub("Bin","",ct1_gr$binID[length(ct1_gr)]))
  rSamp <- sort(rep(sample.int(maxBin,size = 1),length(b)) + b)
  sampBins <- paste0("Bin",rSamp)
  
  binned_sampled_wSE <- left_join(binnedGM %>% filter(binID %in% sampBins),SE_genes, by = "Gene")
  print(dim(seList[[i+1]]))
  toSample <- nrow(seList[[i+1]])
  SE_null[[i+1]] <- binned_sampled_wSE %>% drop_na() %>% ungroup() %>% 
    slice_sample(n = toSample) %>% dplyr::select(SE) %>% mutate(Cluster = i)
}

nullDF <- do.call('rbind',SE_null)
nullDF$Type <- "Null"

fullDF <- rbind(seDF,nullDF)
fullDF$Cluster <- as.factor(fullDF$Cluster)

ggplot(fullDF, aes(x = Cluster, y = SE, color = Type)) +
  geom_boxplot() +ylim(0.75,1)+ 
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("grey60","grey20")) +
  theme(legend.position = "bottom")

ggplot(fullDF, aes(x = Type, y = SE, color = Type)) +
  geom_boxplot() +ylim(0.85,1.1) + 
  geom_signif(comparisons = list(c("Null","Weighted"))) +
  facet_wrap(~Cluster)
