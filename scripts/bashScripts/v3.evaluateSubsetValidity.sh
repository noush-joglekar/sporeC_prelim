#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=eval_subset
#SBATCH --time=04:00:00   # HH/MM/SS
#SBATCH --mem=16G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date` >> evalSubset.log

echo "This is job number:" $SLURM_JOB_ID >> evalSubset.log
echo "Running on node:" `hostname` >> evalSubset.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> evalSubset.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> evalSubset.log

source activate hypergraph_poreC

echo "Trial run for 200 reads"

dirName=$1;

fileName=$(ls $dirName/*.pkl)

echo $dirName;
echo $fileName;

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/getMultiwayExpectedProbsByDist.py \
v2.checkSubsetValidity/$dirName/ Plots/ dfs/ \
v2.checkSubsetValidity/$fileName \
probHash 101 200 25 --plotRef --plotScatter &>> evalSubset.log

echo "Run finished for 200 reads per card and scatterplots generated"
echo "Starting with full dataset now"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/getMultiwayExpectedProbsByDist.py \
v2.checkSubsetValidity/$dirName/ Plots/ dfs/ \
v2.checkSubsetValidity/$fileName \
probHash 101 2000000 25 &>> evalSubset.log


echo "Run finished for all reads per card"

echo "Job finished or was terminated, please check logs"
exit
