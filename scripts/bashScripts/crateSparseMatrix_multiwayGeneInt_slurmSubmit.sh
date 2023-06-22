#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=multiway_SE_sparseMat
#SBATCH --time=4:00:00   # HH/MM/SS
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

module load R/4.0.0

echo "Starting at:" `date` >> multiway_SE_sparseMat.log

echo "This is job number:" $SLURM_JOB_ID >> multiway_SE_sparseMat.log
echo "Running on node:" `hostname` >> multiway_SE_sparseMat.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> multiway_SE_sparseMat.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> multiway_SE_sparseMat.log

Rscript /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/createSparseMatrix_multiwayGeneInt.R

echo "Job finished or was terminated, please check logs"
exit
