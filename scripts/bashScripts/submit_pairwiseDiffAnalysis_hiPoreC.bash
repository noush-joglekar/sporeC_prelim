#!/bin/bash

echo "Submitting for HiPoreC cell lines"
for card in 3 4 5 ; do echo $card; sbatch \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/v1.runDifferentialAnalysis_slurm_HiPoreC.sh \
differentialAnalysis/HiPoreC_cL/ K562 GM12878 $card ; done
