#!/bin/bash

#SBATCH --partition=bigmem   # cluster-specific
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=arrSynPerChr
#SBATCH --time=24:00:00   # HH/MM/SS
#SBATCH --mem=300G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> synergyPerChrom.log

echo "This is job number:" $SLURM_JOB_ID >> synergyPerChrom.log
echo "Running on node:" `hostname` >> synergyPerChrom.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> synergyPerChrom.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> synergyPerChrom.log

module load R/4.1.2

configFile=$1;

chrID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $configFile);
outDir=$2;

inFile='NlaIII_GM12878_output_byChr/NlaIII_GM12878_'$chrID'.gz'

echo "Processing chr"$chrID

## IMPORTANT: Working directory is set in the R script so give paths relative to that!

Rscript \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/runSynergyModel_byChr_Prometh.R \
$inFile $chrID $outDir

echo "Job finished or was terminated, please check logs"
exit
