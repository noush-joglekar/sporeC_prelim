#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=1
#SBATCH --job-name=eval_GM12878
#SBATCH --time=03:00:00   # HH/MM/SS
#SBATCH --mem=48G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> eval_gm12878.log

echo "This is job number:" $SLURM_JOB_ID >> eval_gm12878.log
echo "Running on node:" `hostname` >> eval_gm12878.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> eval_gm12878.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> eval_gm12878.log

echo "Processing "$SLURM_ARRAY_TASK_ID

chr="chr"$SLURM_ARRAY_TASK_ID

source activate hypergraph_poreC

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
v1.evaluateExpectedVersusInteresting_NlaIII_GM12878/ Plots_${chr}/ dfs_${chr}/ \
NlaIII_GM12878_output_byChr/NlaIII_GM12878_${chr}.gz ${chr} \
hyperEdges_GM12878_${chr}.pkl probHash_GM12878_${chr} \
101 3 25 --plotInd --plotRef &>> eval_gm12878.log

echo "Run finished for 3 reads per card and individual reads plotted"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
v1.evaluateExpectedVersusInteresting_NlaIII_GM12878/ Plots_${chr}/ dfs_${chr}/ \
NlaIII_GM12878_output_byChr/NlaIII_GM12878_${chr}.gz ${chr} \
hyperEdges_GM12878_${chr}.pkl probHash_GM12878_${chr} \
101 200 25 --plotScatter &>> eval_gm12878.log

echo "Run finished for 200 reads per card and scatterplots generated"
echo "Starting with full dataset now"

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/promethData_evaluateExpectedVersusInteresting.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v1_poreC_explore/ \
v1.evaluateExpectedVersusInteresting_NlaIII_GM12878/ Plots_${chr}/ dfs_${chr}/ \
NlaIII_GM12878_output_byChr/NlaIII_GM12878_${chr}.gz ${chr} \
hyperEdges_GM12878_${chr}.pkl probHash_GM12878_${chr} \
101 2000000 25 &>> eval_gm12878.log

echo "Run finished for 2 mil reads per card"

echo "Job finished or was terminated, please check logs"
exit
