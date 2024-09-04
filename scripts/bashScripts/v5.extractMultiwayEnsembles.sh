#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --array=1-10000%500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=10k_both
#SBATCH --time=03:00:00   # HH/MM/SS
#SBATCH --mem=30G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

source activate hypergraph_poreC

fileID=$1;

echo "Processing "$SLURM_ARRAY_TASK_ID

## IMPORTANT: Directories hard coded in python script

echo "Cell type 1"
time python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/extractEnsembles_multiwayConstrained.py \
$SLURM_ARRAY_TASK_ID inputData/cell_1_10000_chains_reconstructed_800/ cell1_output/ --prim_cutoff 150 --sec_cutoff 300 \
--num_processes 4 --offDiagDist 1 --cellType A

echo "Cell type 2"
time python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/extractEnsembles_multiwayConstrained.py \
$SLURM_ARRAY_TASK_ID inputData/cell_2_10000_chains_reconstructed_800/ cell2_output/ --prim_cutoff 150 --sec_cutoff 300 \
--num_processes 4 --offDiagDist 1 --cellType B

echo "Job finished or was terminated, please check logs"
exit
