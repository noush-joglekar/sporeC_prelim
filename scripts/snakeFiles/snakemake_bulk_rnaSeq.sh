#!/bin/bash

export PATH=$PATH:/nfs/sw/miniconda3/miniconda3-3.22.0/bin
module load picard/2.8.0

snakemake -d /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/2023_03_06_GM12878_cellularFractionData/ \
    --config ymls=$(pwd)/ymls \
    --configfile=snakeFileConfigs/config_bulk_rnaSeq.yaml \
    --use-conda \
    --conda-frontend conda \
    --profile profile_bulk_rnaSeq \
    --keep-going \
    -s Snakefile_bulk_rnaSeq_new \
    "$@"
