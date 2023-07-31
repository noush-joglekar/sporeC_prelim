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

echo "Starting at:" `date` >> promethPoreC_preprocess.log

echo "This is job number:" $SLURM_JOB_ID >> promethPoreC_preprocess.log
echo "Running on node:" `hostname` >> promethPoreC_preprocess.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> promethPoreC_preprocess.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> promethPoreC_preprocess.log

for i in $(ls NlaIII_GM12878_data/*alignments.csv.gz) ; do batch=$(echo $i | awk '{split($1,b,"_"); print b[6]}') ; echo "Processing batch "$batch; \
sh /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/getGeneNamesFromChromunities_v1.bash $i \
$batch NlaIII_GM12878_output/ ; done

rm fragFile*

echo "Job finished or was terminated, please check logs"
exit
