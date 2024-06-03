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

CT1=$1;
CT2=$2;
card=$3;

chr="chr"$SLURM_ARRAY_TASK_ID
echo "Processing "$SLURM_ARRAY_TASK_ID

echo "About to visualize" $CT1 "versus" $CT2 "for card ="$card " and chrom = "$chr

Rscript /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/plot_ctSpecificConcatemers_hiPoreC.R \
'/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/HiPoreC' \
'differentialAnalysis/HiPoreC_cL/' $CT1 $CT2 $card $chr 'differentialAnalysis/Plots_K562_vs_GM12878/' \
'compartmentData/K562/K562_genomeComp.bigWig' 'compartmentData/GM12878/GM12878_genomeComp.bigWig'
