#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=evalDistCutoff
#SBATCH --time=01:00:00   # HH/MM/SS
#SBATCH --mem=2G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date` >> evalDistCutoff.log

echo "This is job number:" $SLURM_JOB_ID >> evalDistCutoff.log
echo "Running on node:" `hostname` >> evalDistCutoff.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> evalDistCutoff.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> evalDistCutoff.log

source activate hypergraph_poreC

outDir=$1;
jobID=$2;
distCutoff_prim=$3;
distCutoff_sec=$4;
numThreads=$5;
offDiagLim=$6;

echo "Processing "$jobID

## IMPORTANT: Directories hard coded in python script

time python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/makeProjectionMatrices.py \
$outDir $jobID $distCutoff_prim $distCutoff_sec $numThreads --offDiagLim $offDiagLim

echo "Job finished or was terminated, please check logs"
exit
