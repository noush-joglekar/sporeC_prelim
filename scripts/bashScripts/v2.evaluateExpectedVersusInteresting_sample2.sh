#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=eval_S2_10k
#SBATCH --time=06:00:00   # HH/MM/SS
#SBATCH --mem=8G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date` >> eval_10k_sample2.log

echo "This is job number:" $SLURM_JOB_ID >> eval_10k_sample2.log
echo "Running on node:" `hostname` >> eval_10k_sample2.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> eval_10k_sample2.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> eval_10k_sample2.log

source activate hypergraph_poreC

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/getMultiwayExpectedProbsByDist.py \
v2.evaluateExpectedVersusInteresting_sample2/ Plots_sample2_10kChains_500_750_1/ dfs_sample2_10kChains_500_750_1/ \
v2.makeCombinedHypergraphDicts/sample2/hyperEdges_1_500_750_final_chains.pkl \
probHash_10kChains_1_500_750 101 3 25 --plotInd --plotRef &>> eval_10k_sample2.log

echo "Run finished for 3 reads per card and individual reads plotted"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/getMultiwayExpectedProbsByDist.py \
v2.evaluateExpectedVersusInteresting_sample2/ Plots_sample2_10kChains_500_750_1/ dfs_sample2_10kChains_500_750_1/ \
v2.makeCombinedHypergraphDicts/sample2/hyperEdges_1_500_750_final_chains.pkl \
probHash_10kChains_1_500_750 101 200 25 --plotScatter &>> eval_10k_sample2.log

echo "Run finished for 200 reads per card and scatterplots generated"
echo "Starting w/ full dataset"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/getMultiwayExpectedProbsByDist.py \
v2.evaluateExpectedVersusInteresting_sample2/ Plots_sample2_10kChains_500_750_1/ dfs_sample2_10kChains_500_750_1/ \
v2.makeCombinedHypergraphDicts/sample2/hyperEdges_1_500_750_final_chains.pkl \
probHash_10kChains_1_500_750 101 2000000 25 &>> eval_10k_sample2.log

echo "Run finished for 2 mil reads per card"

echo "Job finished or was terminated, please check logs"
exit
