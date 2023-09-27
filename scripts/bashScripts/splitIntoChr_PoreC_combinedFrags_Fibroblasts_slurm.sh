#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fib_splitChr
#SBATCH --time=00:05:00   # HH/MM/SS
#SBATCH --mem=500M   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

echo "Starting at:" `date` >> fib_splitChr.log.log

echo "This is job number:" $SLURM_JOB_ID >> fib_splitChr.log.log
echo "Running on node:" `hostname` >> fib_splitChr.log.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> fib_splitChr.log.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> fib_splitChr.log.log

chrID=$1;
file=$2;
prefix=$3;

echo "Processing chr"$chrID

awk -v chr="chr"${chrID} -v file=$file 'BEGIN{c=1; comm="zcat "file; \
while(comm|getline) {if($1==chr) {if(!id[$4]) {id[$4]=c; c=c+1;} $15=id[$4];
OFS="\t"; print }} }' | gzip -c > ${prefix}_chr${chrID}.gz


echo "Job finished or was terminated, please check logs"
exit
