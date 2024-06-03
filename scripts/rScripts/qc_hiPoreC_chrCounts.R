## QC on HiPoreC cell lines
library(ggplot2)

setwd('../HiPoreC/')
chrNames <- paste("chr",1:22,sep = "")

k562 <- read.table('DpnII_K562_output/K562_chrCounts',sep = "\t")
gm <- read.table('DpnII_GM12878_output/GM12878_chrCounts',sep = "\t")
colnames(k562) <- colnames(gm) <- c("chrName","counts")
k562 <- k562 %>% filter(chrName %in% chrNames) %>% mutate(CL = "K562")
gm <- gm %>% filter(chrName %in% chrNames) %>% mutate(CL = "GM12878")

fullDF <- rbind(k562,gm)
fullDF$chrName <- factor(fullDF$chrName, levels = chrNames)


ggplot(fullDF, aes(x = chrName, y = counts, fill = CL)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~CL) +
  theme(axis.text.x = element_text(angle = 90))
