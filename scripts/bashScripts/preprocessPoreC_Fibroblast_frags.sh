#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=procFibs
#SBATCH --time=02:00:00   # HH/MM/SS
#SBATCH --mem=4G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL


echo "Starting at:" `date` >> fibroPoreC_preprocess.log

echo "This is job number:" $SLURM_JOB_ID >> fibroPoreC_preprocess.log
echo "Running on node:" `hostname` >> fibroPoreC_preprocess.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> fibroPoreC_preprocess.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> fibroPoreC_preprocess.log

source activate hypergraph_poreC

sample=$1;

python /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/parquet2csv.py \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/PoreC_Fibroblasts/ \
${sample}/ ${sample}_csv/ chr_alias.tbl

for i in $(ls ${sample}_csv/*.csv.gz) ; do batch=$(echo $i | awk '{split($1,b,"_"); print b[3]}') ; echo "Processing batch "$batch; \
sh /gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/bashScripts/getGeneNamesFromChromunities_v1.bash $i \
$batch ${sample}_fragsOutput/ ; done

rm fragFile*

echo "Job finished or was terminated, please check logs"
exit
