library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(parallel)

## Wrapping in bash for slurm submission
## Associated script: /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/crateSparseMatrix_multiwayGeneInt_slurmSubmit.sh

scriptDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/'
workingDir <- '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_poreC_SE_interaction/'
setwd(workingDir)

readChromDF <- fread('../v0_poreC_explore/readChromunities_wClosestGene.gz')
colnames(readChromDF) <- c("ReadID","Cardinality","gene_ID","geneName","distToGene")

## Create a sparse matrix
geneList <- unique(readChromDF$geneName)
nGenes <- length(geneList)

print(paste(nGenes,"total genes"))

uniqueReads <- unique(readChromDF$ReadID)
numBins <- ceiling(length(uniqueReads)/5000) ##10000
chunkGenes <- split(uniqueReads, 1:numBins)
print(paste0(numBins," bins generated"))

getGenePairs <- function(readID,DF,sparseM){
  subset <- DF %>% filter(ReadID == readID)
  geneIndices <- match(subset$geneName, geneList)
  geneCombs <- combn(geneIndices,2) ## can modify to make it only unique genes in the read chromunity
  sparseM[geneCombs[1, ], geneCombs[2, ]] <- sparseM[geneCombs[1, ], geneCombs[2, ]] + 1
  sparseM[geneCombs[2, ], geneCombs[1, ]] <- sparseM[geneCombs[2, ], geneCombs[1, ]] + 1
  return(sparseM)
}

f0 = Sys.time()
for(c in 1:numBins){
  print(c)
  chunk <- chunkGenes[[c]]
  DF <- readChromDF %>% filter(ReadID %in% chunk)
  sparseM <- Matrix(0,nGenes,nGenes,sparse = TRUE)
  A = mclapply(chunk,function(r) getGenePairs(r,DF,sparseM),mc.cores = 24)
  sparseM <- Reduce('+', A)
  if(file.exists('tmpMatrix.mtx')){
    oldSM <- readMM('tmpMatrix.mtx')
    newSM <- oldSM + sparseM
    writeMM(newSM,'tmpMatrix.mtx')
    rm(A,newSM, oldSM,sparseM)
  } else {
    writeMM(sparseM,'tmpMatrix.mtx')
    rm(A,sparseM)
  }
}
f1 = Sys.time()
f1-f0

print("Final wrap-up in progress")
fullSM <- readMM('tmpMatrix.mtx')
colnames(fullSM) <- rownames(fullSM) <- geneList
writeMM(fullSM,'allGenesMatrix.mtx')
file.remove('tmpMatrix.mtx')

print("All done!")
