#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --array=1-10000%500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=arrPM
#SBATCH --time=00:45:00   # HH/MM/SS
#SBATCH --mem=1500M   # memory requested, units available: K,M,G,T
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
chains_10k_500_projectionMtxOutput/ $SLURM_ARRAY_TASK_ID 600 750 3 --offDiagLim 3

echo "Job finished or was terminated, please check logs"
exit
