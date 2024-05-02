#!/bin/bash

echo "Submitting for HCC versus HG"
for card in 3 4 5 ; do echo $card; sbatch \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/v1.runDifferentialAnalysis_slurm.sh \
differentialAnalysis/hcc_vs_hg/ HCC1954 HG002 $card ; done

echo "Submitting for HCC versus GM"
for card in 3 4 5 ; do echo $card; sbatch \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/v1.runDifferentialAnalysis_slurm.sh \
differentialAnalysis/gm_vs_hcc/ HCC1954 GM12878 $card ; done
