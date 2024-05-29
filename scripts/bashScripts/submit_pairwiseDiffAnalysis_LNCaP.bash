#!/bin/bash

echo "Submitting for LNCaP"
for card in 3 4 5 ; do echo $card; sbatch \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/v1.runDifferentialAnalysis_slurm_LNCaP.sh \
differentialAnalysis/LNCaP/ LNCaP-Vehicle LNCaP-DHT $card ; done
