## QC on cell lines
library(ggplot2)

chrNames <- paste("chr",1:22,sep = "")

## PRE PROCESSING ----
hcc_pre1 <- read.table('NlaIII_HCC1954_data/3c36a_chrCounts',sep = "\t")
hcc_pre2 <- read.table('NlaIII_HCC1954_data/6de06_chrCounts',sep = "\t")
hcc_pre <- rbind(hcc_pre1,hcc_pre2) %>% group_by(V1) %>% mutate(V2 = sum(V2)) %>% distinct()
colnames(hcc_pre) <- c("chrName","counts")
hcc_pre <- hcc_pre %>% filter(chrName %in% chrNames) %>% mutate(CL = "HCC1954", Type = "Pre")

hg_pre1 <- read.table('NlaIII_HG002_data/bdb7b_chrCounts',sep = "\t")
hg_pre2 <- read.table('NlaIII_HG002_data/1c575_chrCounts',sep = "\t")
hg_pre <- rbind(hg_pre1,hg_pre2) %>% group_by(V1) %>% mutate(V2 = sum(V2)) %>% distinct()
colnames(hg_pre) <- c("chrName","counts")
hg_pre <- hg_pre %>% filter(chrName %in% chrNames) %>% mutate(CL = "HG002", Type = "Pre")

lnVeh_pre1 <- read.table('NlaIII_LNCaP-Vehicle_data/PAG48791_chrCounts',sep = "\t")
lnVeh_pre2 <- read.table('NlaIII_LNCaP-Vehicle_data/PAG49193_chrCounts',sep = "\t")
lnVeh_pre <- rbind(lnVeh_pre1,lnVeh_pre2) %>% group_by(V1) %>% mutate(V2 = sum(V2)) %>% distinct()
colnames(lnVeh_pre) <- c("chrName","counts")
lnVeh_pre <- lnVeh_pre %>% filter(chrName %in% chrNames) %>% mutate(CL = "LnCap_Vehicle", Type = "Pre")

lnDHT_pre1 <- read.table('NlaIII_LNCaP-DHT_data/PAG49593_chrCounts',sep = "\t")
lnDHT_pre2 <- read.table('NlaIII_LNCaP-DHT_data/PAG49696_chrCounts',sep = "\t")
lnDHT_pre <- rbind(lnDHT_pre1,lnDHT_pre2) %>% group_by(V1) %>% mutate(V2 = sum(V2)) %>% distinct()
colnames(lnDHT_pre) <- c("chrName","counts")
lnDHT_pre <- lnDHT_pre %>% filter(chrName %in% chrNames) %>% mutate(CL = "LnCap_DHT", Type = "Pre")

## POST PROCESSING ----

hcc_post <- read.table('postProcessingChrCounts/NlaIII_HCC1954_output_chrCounts',sep = "\t")
colnames(hcc_post) <- c("chrName","counts")
hcc_post <- hcc_post %>% filter(chrName %in% chrNames) %>% mutate(CL = "HCC1954", Type = "Post")

hg_post <- read.table('postProcessingChrCounts/NlaIII_HG002_output_chrCounts',sep = "\t")
colnames(hg_post) <- c("chrName","counts")
hg_post <- hg_post %>% filter(chrName %in% chrNames) %>% mutate(CL = "HG002", Type = "Post")

lnVeh_post <- read.table('postProcessingChrCounts/NlaIII_LNCaP-Vehicle_output_chrCounts',sep = "\t")
colnames(lnVeh_post) <- c("chrName","counts")
lnVeh_post <- lnVeh_post %>% filter(chrName %in% chrNames) %>% mutate(CL = "LnCap_Vehicle", Type = "Post")

lnDHT_post <- read.table('postProcessingChrCounts/NlaIII_LNCaP-DHT_output_chrCounts',sep = "\t")
colnames(lnDHT_post) <- c("chrName","counts")
lnDHT_post <- lnDHT_post %>% filter(chrName %in% chrNames) %>% mutate(CL = "LnCap_DHT", Type = "Post")

fullDF <- rbind(hcc_pre,hg_pre,lnVeh_pre,lnDHT_pre,hcc_post, hg_post, lnVeh_post,lnDHT_post)
fullDF$chrName <- factor(fullDF$chrName, levels = chrNames)

fullDF$Type <- factor(fullDF$Type, levels <- c("Pre","Post"))
pdf('postProcessingChrCounts/QC_allCellLines_preAndPost_corrected.pdf',10,10,useDingbats = F)
ggplot(fullDF, aes(x = chrName, y = counts, fill = CL)) + 
  geom_bar(stat = "identity") +
  facet_grid(Type~CL) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
