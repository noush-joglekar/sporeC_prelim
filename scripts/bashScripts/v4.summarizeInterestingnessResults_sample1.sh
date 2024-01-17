#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=summ_s1
#SBATCH --time=04:00:00   # HH/MM/SS
#SBATCH --mem=12G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> summarizeResults_sample1.log

echo "This is job number:" $SLURM_JOB_ID >> summarizeResults_sample1.log
echo "Running on node:" `hostname` >> summarizeResults_sample1.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> summarizeResults_sample1.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> summarizeResults_sample1.log

source activate hypergraph_poreC

echo "Starting job ...."

sampleSheet=$1;

while read line; do \
dname=$(echo $line | awk '{print $4}'); \
echo $line > sample1/$dname/tmpFile1; \
time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/checkInterestingnessStatusOfSubsets.py \
sample1/$dname/ summaryFile_sample1.tab --readOfInterest tmpFile1 ; \
rm sample1/$dname/tmpFile1; \
done < $sampleSheet

echo "Job finished or was terminated, please check logs"
exit
