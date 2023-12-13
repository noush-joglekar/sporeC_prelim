#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=chunkPklDict_s1
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=13G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> combineAllChains_dict_sample1.log

echo "This is job number:" $SLURM_JOB_ID >> combineAllChains_dict_sample1.log
echo "Running on node:" `hostname` >> combineAllChains_dict_sample1.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> combineAllChains_dict_sample1.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> combineAllChains_dict_sample1.log

source activate hypergraph_poreC

## IMPORTANT: Directories hard coded in python script

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/v1.combineAllChains.py \
10000 500 v2.processChainsOutput_10k_500_sample1/ v2.makeCombinedHypergraphDicts/sample1/ \
--prim_cutoff 550 --sec_cutoff 750 --offDiagDist 3 --num_cores 8

echo "Job finished or was terminated, please check logs"
exit
