#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=comb
#SBATCH --time=18:00:00   # HH/MM/SS
#SBATCH --mem=12G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date`

echo "This is job number:" $SLURM_JOB_ID
echo "Running on node:" `hostname`
echo "Running on cluster:" $SLURM_CLUSTER_NAME
echo "This job was assigned the temporary (local) directory:" $TMPDIR

source activate hypergraph_poreC

## IMPORTANT: Directories hard coded in python script

echo "Running cell 1"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/v1.combineAllChains_multiConstrained.py \
1000 500 cell1_output/ comb_hpDicts/cell1/ \
--prim_cutoff 150 --sec_cutoff 300 --offDiagDist 1 --num_cores 8

echo "Running cell 2"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/v1.combineAllChains_multiConstrained.py \
1000 500 cell2_output/ comb_hpDicts/cell2/ \
--prim_cutoff 150 --sec_cutoff 300 --offDiagDist 1 --num_cores 8
