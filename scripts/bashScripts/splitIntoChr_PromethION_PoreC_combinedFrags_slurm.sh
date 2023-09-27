#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=splitChr_prom
#SBATCH --time=00:30:00   # HH/MM/SS
#SBATCH --mem=5G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date` >> splitChr_prom.log

echo "This is job number:" $SLURM_JOB_ID >> splitChr_prom.log
echo "Running on node:" `hostname` >> splitChr_prom.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> splitChr_prom.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> splitChr_prom.log

chrID=$1;
file=$2;
prefix=$3;

echo "Processing chr"$chrID

awk -v chr="chr"${chrID} -v file=$file 'BEGIN{c=1; comm="zcat "file; \
while(comm|getline) {if($1==chr) {if(!id[$4]) {id[$4]=c; c=c+1;} $15=id[$4];
OFS="\t"; print }} }' | gzip -c > ${prefix}_chr${chrID}.gz

echo "Job finished or was terminated, please check logs"
exit
