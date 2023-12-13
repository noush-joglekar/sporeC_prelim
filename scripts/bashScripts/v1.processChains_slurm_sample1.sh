#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --array=1-10000%250
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=10k_s1
#SBATCH --time=03:00:00   # HH/MM/SS
#SBATCH --mem=6G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date` >> 10k_sample1_550_750_3.log

echo "This is job number:" $SLURM_JOB_ID >> 10k_sample1_550_750_3.log
echo "Running on node:" `hostname` >> 10k_sample1_550_750_3.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> 10k_sample1_550_750_3.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> 10k_sample1_550_750_3.log

source activate hypergraph_poreC

echo "Processing "$SLURM_ARRAY_TASK_ID

## IMPORTANT: Directories hard coded in python script

time python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/v1.processChains.py \
chains_500_10000_1500_1681171613/ v2.processChainsOutput_10k_500_sample1/ $SLURM_ARRAY_TASK_ID 550 750 4 --offDiagLim 3

echo "Job finished or was terminated, please check logs"
exit
