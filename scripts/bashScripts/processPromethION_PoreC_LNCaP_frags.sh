#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=promethPoreC
#SBATCH --time=8:00:00   # HH/MM/SS
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

## Had to make modifications for LNCaP

inputDir=$1;
outputDir=$2;
logFile=$3;

echo "Starting at:" `date` >> $logFile

echo "This is job number:" $SLURM_JOB_ID >> $logFile
echo "Running on node:" `hostname` >> $logFile
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> $logFile
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> $logFile

for i in $(ls $inputDir/*full_pore_c.csv.gz) ; do batch=$(echo $i | awk '{split($1,b,"_"); print b[7]}') ; echo "Processing batch "$batch; \
sh /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/getGeneNamesFromChromunities_LNCaP.bash $i \
$batch $outputDir ; done


echo "Job finished or was terminated, please check logs"
exit
