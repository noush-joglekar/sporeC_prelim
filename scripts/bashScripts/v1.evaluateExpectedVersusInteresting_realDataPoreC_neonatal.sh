#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --job-name=eval_neonatal
#SBATCH --time=1:00:00   # HH/MM/SS
#SBATCH --mem=8G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> eval_neonatal.log

echo "This is job number:" $SLURM_JOB_ID >> eval_neonatal.log
echo "Running on node:" `hostname` >> eval_neonatal.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> eval_neonatal.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> eval_neonatal.log

echo "Processing "$SLURM_ARRAY_TASK_ID

chr="chr"$SLURM_ARRAY_TASK_ID

source activate hypergraph_poreC

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/PoreC_Fibroblasts/ \
v1.evaluateExpectedVersusInteresting_neonatal/ Plots_${chr}/ dfs_${chr}/ \
neonatal_fragsOutput_byChr/neonatal_fibroblasts_${chr}.gz ${chr} \
hyperEdges_neonatal_${chr}.pkl probHash_neonatal_${chr} \
101 3 25 --plotInd --plotRef &>> eval_neonatal.log

echo "Run finished for 3 reads per card and individual reads plotted"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/PoreC_Fibroblasts/ \
v1.evaluateExpectedVersusInteresting_neonatal/ Plots_${chr}/ dfs_${chr}/ \
neonatal_fragsOutput_byChr/neonatal_fibroblasts_${chr}.gz ${chr} \
hyperEdges_neonatal_${chr}.pkl probHash_neonatal_${chr} \
101 200 25 --plotScatter &>> eval_neonatal.log

echo "Run finished for 200 reads per card and scatterplots generated"
echo "Starting with full dataset now"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/PoreC_Fibroblasts/ \
v1.evaluateExpectedVersusInteresting_neonatal/ Plots_${chr}/ dfs_${chr}/ \
neonatal_fragsOutput_byChr/neonatal_fibroblasts_${chr}.gz ${chr} \
hyperEdges_neonatal_${chr}.pkl probHash_neonatal_${chr} \
101 2000000 25 &>> eval_neonatal.log

echo "Run finished for 2 mil reads per card"

echo "Job finished or was terminated, please check logs"
exit
