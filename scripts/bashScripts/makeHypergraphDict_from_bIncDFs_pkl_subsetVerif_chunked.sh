#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=subsetVerif_iter3
#SBATCH --time=10:00:00   # HH/MM/SS
#SBATCH --mem=12G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=ALL

echo "Starting at:" `date` >> hpDict.log

echo "This is job number:" $SLURM_JOB_ID >> hpDict.log
echo "Running on node:" `hostname` >> hpDict.log
echo "Running on cluster:" $SLURM_CLUSTER_NAME >> hpDict.log
echo "This job was assigned the temporary (local) directory:" $TMPDIR >> hpDict.log

source activate hypergraph_poreC

numFiles=$1;
offDiagDist=$2;
inputDir=$3;
outputDir=$4;
chunkSize=$5;
readOfInterest=$6
sampleFract=$7;
## IMPORTANT: Directories hard coded in python script

time python \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/makeHypergraphDict_fromChainIncDFs_pkl_subsetVerif_chunked.py \
$numFiles $offDiagDist $inputDir $outputDir $chunkSize $readOfInterest $sampleFract

echo "Job finished or was terminated, please check logs"
exit
