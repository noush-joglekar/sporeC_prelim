#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --job-name=diffAn_hpc
#SBATCH --time=02:00:00   # HH/MM/SS
#SBATCH --mem=80G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date`

echo "This is job number:" $SLURM_JOB_ID
echo "Running on node:" `hostname`
echo "Running on cluster:" $SLURM_CLUSTER_NAME
echo "This job was assigned the temporary (local) directory:" $TMPDIR

source activate hypergraph_poreC

outDir=$1;
CT1=$2;
CT2=$3;
card=$4;

chr="chr"$SLURM_ARRAY_TASK_ID
echo "Processing "$SLURM_ARRAY_TASK_ID

dirName="/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/HiPoreC/"${outDir}

if [ ! -d $dirName ]; then
    echo "Creating dir"
    mkdir $dirName
fi

echo "About to compare" $CT1 "versus" $CT2 "for card ="$card " and chrom = "$chr

python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/v4.differentialAnalysis.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/HiPoreC/ \
${outDir}/ \
projMatPlots_cellLines_hpc/matrices/ $chr $card $CT1 $CT2 --re DpnII
