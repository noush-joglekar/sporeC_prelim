#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --array=1-1000%100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=arrPM
#SBATCH --time=01:00:00   # HH/MM/SS
#SBATCH --mem=4G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> projMat.log

echo "This is job number:" $SLURM_JOB_ID >> projMat.log
echo "Running on node:" `hostname` >> projMat.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> projMat.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> projMat.log

source activate hypergraph_poreC

echo "Processing "$SLURM_ARRAY_TASK_ID

## IMPORTANT: Directories hard coded in python script

time python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/makeProjectionMatrices.py \
chains_10k_500_projectionMtxOutput_v2/ $SLURM_ARRAY_TASK_ID 600 750 4 --offDiagLim 1

echo "Job finished or was terminated, please check logs"
exit
