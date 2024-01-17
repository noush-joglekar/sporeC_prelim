#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=vc_s1
#SBATCH --time=08:00:00   # HH/MM/SS
#SBATCH --mem=12G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> validityCheck_sample1.log

echo "This is job number:" $SLURM_JOB_ID >> validityCheck_sample1.log
echo "Running on node:" `hostname` >> validityCheck_sample1.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> validityCheck_sample1.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> validityCheck_sample1.log

source activate hypergraph_poreC

echo "Starting job ...."

sampleSheet=$1;

while read line; do \
dname=$(echo $line | awk '{print $4}'); \
mkdir -p sample1/$dname/ ; \
echo $line > sample1/$dname/tmpFile1; \
time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/extractSubsetChains_forInterestingReads.py \
10000 1000 v2.processChainsOutput_10k_500_sample1/ v2.checkSubsetValidity/sample1/$dname/ \
--prim_cutoff 500 --sec_cutoff 750 --offDiagDist 1 --num_cores 8 --readOfInterest tmpFile1 ; \
rm sample1/$dname/tmpFile1; \
done < $sampleSheet

echo "Job finished or was terminated, please check logs"
exit
