#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --job-name=GM_3_50KB
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

cellLine=$1;
minReads=$2;
## Only works with GM12878 and minimum reads = 2,3. Bin size is hardcoded at 50KB

echo "Starting at:" `date` >> eval_${cellLine}.log

echo "This is job number:" $SLURM_JOB_ID >> eval_${cellLine}.log
echo "Running on node:" `hostname` >> eval_${cellLine}.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> eval_${cellLine}.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> eval_${cellLine}.log

echo "Processing "$SLURM_ARRAY_TASK_ID

chr="chr"$SLURM_ARRAY_TASK_ID

source activate hypergraph_poreC

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
v1.evaluateExpectedVersusInteresting_NlaIII_50KB_${cellLine}_${minReads}/ Plots_${chr}/ dfs_${chr}/ \
NlaIII_${cellLine}_output_byChr/NlaIII_${cellLine}_${chr}.gz ${chr} \
hyperEdges_${cellLine}_${chr}.pkl probHash_${cellLine}_${chr} \
101 3 25 --plotInd --plotRef

echo "Run finished for 3 reads per card and individual reads plotted"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
v1.evaluateExpectedVersusInteresting_NlaIII_50KB_${cellLine}_${minReads}/ Plots_${chr}/ dfs_${chr}/ \
NlaIII_${cellLine}_output_byChr/NlaIII_${cellLine}_${chr}.gz ${chr} \
hyperEdges_${cellLine}_${chr}.pkl probHash_${cellLine}_${chr} \
101 200 25 --plotScatter

echo "Run finished for 200 reads per card and scatterplots generated"
echo "Starting with full dataset now"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
v1.evaluateExpectedVersusInteresting_NlaIII_50KB_${cellLine}_${minReads}/ Plots_${chr}/ dfs_${chr}/ \
NlaIII_${cellLine}_output_byChr/NlaIII_${cellLine}_${chr}.gz ${chr} \
hyperEdges_${cellLine}_${chr}.pkl probHash_${cellLine}_${chr} \
101 2000000 25

echo "Run finished for 2 mil reads per card"

echo "Job finished or was terminated, please check logs"
exit
