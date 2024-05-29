#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-24
#SBATCH --cpus-per-task=1
#SBATCH --job-name=hpc_75
#SBATCH --time=4:00:00   # HH/MM/SS
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

## Works with GM12878, HG002, HCC1954 with the percentage-based min-read cutoff. Bin size is 1MB. Cutoff is 75% (25)
## Should also work with LNCaP
## Running for HiPoreC. K

cellLine=$1;
# minReads=$2;

echo "Starting at:" `date` >> eval_${cellLine}.log

echo "This is job number:" $SLURM_JOB_ID >> eval_${cellLine}.log
echo "Running on node:" `hostname` >> eval_${cellLine}.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> eval_${cellLine}.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> eval_${cellLine}.log

if [[ "$SLURM_ARRAY_TASK_ID" == 23 ]]; then
chr="chrX"
echo "Processing chrX"
elif [[ "$SLURM_ARRAY_TASK_ID" == 24 ]]; then
chr="chrY"
echo "Processing chrY"
else
chr="chr"$SLURM_ARRAY_TASK_ID
echo "Processing "$SLURM_ARRAY_TASK_ID
fi

source activate hypergraph_poreC

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/HiPoreC/ \
v1.evaluateExpectedVersusInteresting_DpnII_${cellLine}/ Plots_${chr}/ dfs_${chr}/ \
DpnII_${cellLine}_output_byChr/DpnII_${cellLine}_${chr}.gz ${chr} \
hyperEdges_${cellLine}_${chr}.pkl probHash_${cellLine}_${chr} \
101 3 25 --plotInd --plotRef

echo "Run finished for 3 reads per card and individual reads plotted"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/HiPoreC/ \
v1.evaluateExpectedVersusInteresting_DpnII_${cellLine}/ Plots_${chr}/ dfs_${chr}/ \
DpnII_${cellLine}_output_byChr/DpnII_${cellLine}_${chr}.gz ${chr} \
hyperEdges_${cellLine}_${chr}.pkl probHash_${cellLine}_${chr} \
101 200 25 --plotScatter

echo "Run finished for 200 reads per card and scatterplots generated"
echo "Starting with full dataset now"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/HiPoreC/ \
v1.evaluateExpectedVersusInteresting_DpnII_${cellLine}/ Plots_${chr}/ dfs_${chr}/ \
DpnII_${cellLine}_output_byChr/DpnII_${cellLine}_${chr}.gz ${chr} \
hyperEdges_${cellLine}_${chr}.pkl probHash_${cellLine}_${chr} \
101 3000000 25

echo "Run finished for 3 mil reads per card"

echo "Job finished or was terminated, please check logs"
exit
