#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=prometh_splitChr
#SBATCH --time=4:00:00   # HH/MM/SS
#SBATCH --mem=16G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> prometh_splitChr_preprocess.log

echo "This is job number:" $SLURM_JOB_ID >> prometh_splitChr_preprocess.log
echo "Running on node:" `hostname` >> prometh_splitChr_preprocess.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> prometh_splitChr_preprocess.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> prometh_splitChr_preprocess.log

chrID=$1;
file=$2;

echo "Processing chr"$chrID

awk -v chr=$chrID -v file=$file 'BEGIN{c=1; comm="zcat "file; \
while(comm|getline) {if($1~chr) {if(!id[$4]) {id[$4]=c; c=c+1;} $15=id[$4];
OFS="\t"; print }} }' | gzip -c > NlaIII_GM12878_chr${chrID}.gz

echo "Job finished or was terminated, please check logs"
exit
