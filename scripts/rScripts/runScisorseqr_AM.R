#!/nfs/sw/modules/R/4.1.2

# By Anoushka Joglekar
# 05.2024

## IN BASH
#source activate hypergraph_poreC
#module load minimap2/2.17
#module load samtools
#module load bedtools

library(scisorseqr)

print("++++++++ Step 1: Getting barcodes and filtering those reads out")
#GetBarcodes(fqFolder='/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2024_05_29_mskiLab_MelanomaIsoSeq/inputData/fastqs/Acral/',
#       BCClustAssignFile='/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2024_05_29_mskiLab_MelanomaIsoSeq/inputData/barcode2cellType',
#       chemistry='v3', filterReads=FALSE, numProcesses=4,concatenate=FALSE)

print("++++++++ Step 2: Aligning with minimap2")
MMalign(fqFolder='/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2024_05_29_mskiLab_MelanomaIsoSeq/inputData/fastqs/Acral/',
        '/nfs/sw/minimap2/minimap2-2.17/minimap2',
        refGenome='/gpfs/commons/groups/gursoy_lab/ajoglekar/Support/References/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly_chr.fa',
        numThreads=4)

print("++++++++ Step 3: Map and filter function")
MapAndFilter(numThreads=4, filterFullLength=TRUE,
        polyABed='/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2024_05_29_mskiLab_MelanomaIsoSeq/inputData/cagePolyA/atlas.clusters.2.0.GRCh38.96_chrNames.bed.gz',
        cageBed='/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2024_05_29_mskiLab_MelanomaIsoSeq/inputData/cagePolyA/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz',
        annoGZ='/gpfs/commons/groups/gursoy_lab/ajoglekar/Support/References/Human/GRCh38/GENCODE/gencode.v42.annotation.gtf.gz',
        seqDir='/gpfs/commons/groups/gursoy_lab/ajoglekar/Support/References/Human/GRCh38/chroms', genomeVersion='hg38')

print("++++++++ Step 4: Getting All-Info files")
InfoPerLongRead(barcodeOutputFile='OutputFiltered/FilteredDeconvBC_AllFiles.csv',
        mapAndFilterOut='LRProcessingOutput/', minTimesIsoObserve=3, rmTmpFolder=FALSE)


#print("++++++++ Step 5: Count isos")
#IsoQuant('LongReadInfo/AllInfo_IncompleteReads.gz')
