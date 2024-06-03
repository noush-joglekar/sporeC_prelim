#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --job-name=diffVis
#SBATCH --time=00:20:00   # HH/MM/SS
#SBATCH --mem=80G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date`

echo "This is job number:" $SLURM_JOB_ID
echo "Running on node:" `hostname`
echo "Running on cluster:" $SLURM_CLUSTER_NAME
echo "This job was assigned the temporary (local) directory:" $TMPDIR

module load R/4.1.2

card=$1;

chr="chr"$SLURM_ARRAY_TASK_ID
echo "Processing "$SLURM_ARRAY_TASK_ID


Rscript /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/plot_ctSpecificConcatemers_cellLines.R \
'/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/' \
'differentialAnalysis/hcc_vs_hg/' HCC1954 HG002 $card $chr 'differentialAnalysis/Plots_HCC_vs_HG/'

Rscript /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/plot_ctSpecificConcatemers_cellLines.R \
'/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/' \
'differentialAnalysis/gm_vs_hcc/' GM12878 HCC1954 $card $chr 'differentialAnalysis/Plots_GM_vs_HCC/'
