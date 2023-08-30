#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=evalDistCutoff
#SBATCH --time=10:00:00   # HH/MM/SS
#SBATCH --mem=4G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> evalDistCutoff.log

echo "This is job number:" $SLURM_JOB_ID >> evalDistCutoff.log
echo "Running on node:" `hostname` >> evalDistCutoff.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> evalDistCutoff.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> evalDistCutoff.log

source activate hypergraph_poreC

jobID=$1;
distCutoff=$2;
numThreads=$3;
offDiagLim=$4;

echo "Processing "$jobID

## IMPORTANT: Directories hard coded in python script

time python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/makeProjectionMatrices.py \
$jobID $distCutoff $numThreads --offDiagLim $offDiagLim

echo "Job finished or was terminated, please check logs"
exit
