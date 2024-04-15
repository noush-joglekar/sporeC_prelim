#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --job-name=plotProj
#SBATCH --time=02:00:00   # HH/MM/SS
#SBATCH --mem=32G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date`

echo "This is job number:" $SLURM_JOB_ID
echo "Running on node:" `hostname`
echo "Running on cluster:" $SLURM_CLUSTER_NAME
echo "This job was assigned the temporary (local) directory:" $TMPDIR

source activate hypergraph_poreC

echo "Starting job ...."

cellLine=$1;
inputDir=$2;
outDir=$3;
empCutoff=$4;
chrom="chr"$SLURM_ARRAY_TASK_ID;

python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/v3.analysisPerChrom.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
${inputDir} ${outDir}/ $chrom $cellLine $empCutoff 3,4,5,6

echo "Job finished or was terminated, please check logs"
exit
