#!/bin/R

## By Anoushka Joglekar 03/2024
## Trying to figure out why the LNCaP data has such a weird distribution of chrom

library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

setwd('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/')

vehicleDir <- 'NlaIII_LNCaP-Vehicle_data/'
dhtDir <- 'NlaIII_LNCaP-DHT_data/'

vFiles <- grep("chrCounts",list.files(vehicleDir),value = T)

vList <- list()
for (file in vFiles){
    name <- unlist(strsplit(file,"_"))[1]
    print(name)
    x <- read.table(file.path(vehicleDir,file))
    colnames(x) <- c("chr","readCounts")
    x$readCounts <- as.numeric(x$readCounts)
    x$Sample <- name
    vList[[name]] <- x
}

vDF <- as.data.frame(do.call('rbind',vList)) %>% remove_rownames()
vDF$Type <- "Vehicle"

dFiles <- grep("chrCounts",list.files(dhtDir),value = T)
dList <- list()
for (file in dFiles){
    name <- unlist(strsplit(file,"_"))[1]
    print(name)
    x <- read.table(file.path(dhtDir,file))
    colnames(x) <- c("chr","readCounts")
    x$readCounts <- as.numeric(x$readCounts)
    x$Sample <- name
    dList[[name]] <- x
}

dDF <- as.data.frame(do.call('rbind',dList)) %>% remove_rownames()
dDF$Type <- "DHT"

fullDF <- rbind(vDF,dDF)
fullDF$chr <- factor(fullDF$chr, levels = paste0("chr",c(1:22,"X","Y","M")))

ggplot(fullDF,aes(x = chr, y = readCounts, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Type) +
    scale_fill_brewer(palette = "Paired") +
    theme_classic(base_size = 14)
