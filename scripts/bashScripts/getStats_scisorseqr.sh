#! /bin/bash
# By Anoushka Joglekar, Tilgner lab
# ca 2022

module load samtools

mb=$(samtools view LRProcessingOutput/mapping.bestperRead.bam | wc -l)
csmm=$(zcat LRProcessingOutput/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz | wc -l)
fl_csmm=$(zcat LRProcessingOutput/CagePolyA.complete.mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz | wc -l)
csmm_b=$(zcat LongReadInfo/AllInfo_IncompleteReads.gz | wc -l)
fl_csmm_b=$(zcat LongReadInfo/AllInfo.gz | wc -l)

echo "Mapped and barcoded= "$mb
echo "CSMM= "$csmm
echo "Full-length CSMM= "$fl_csmm
echo "CSMM barcoded= "$csmm_b
echo "Full-length CSMM barcoded= "$fl_csmm_b
