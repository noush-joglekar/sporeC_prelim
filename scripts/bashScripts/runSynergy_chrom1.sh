#!/bin/bash

#SBATCH --partition=bigmem   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --job-name=synChr2
#SBATCH --time=12:00:00   # HH/MM/SS
#SBATCH --mem=256G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> synergy_chr2.log

echo "This is job number:" $SLURM_JOB_ID >> synergy_chr2.log
echo "Running on node:" `hostname` >> synergy_chr2.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> synergy_chr2.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> synergy_chr2.log

module load R/4.1.2

chrID=$1;
outDir=$2;

inFile='NlaIII_GM12878_output_byChr/NlaIII_GM12878_'$chrID'.gz'

echo "Processing chr"$chrID

## IMPORTANT: Working directory is set in the R script so give paths relative to that!

Rscript \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/rScripts/runSynergyModel_byChr_Prometh.R \
$inFile $chrID $outDir

echo "Job finished or was terminated, please check logs"
exit
