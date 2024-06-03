#!/bin/bash

#SBATCH --partition=pe2   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --job-name=b2f
#SBATCH --time=2:00:00   # HH/MM/SS
#SBATCH --mem=12G   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=ajoglekar@nygenome.org
#SBATCH --mail-type=FAIL

## Takes in file names as input
## Takes in output fastq filename as second argument

module load samtools

inputFile=$1;
outName=$2;

samtools fastq ${inputFile} -@ 2 --reference \
/gpfs/commons/groups/gursoy_lab/ajoglekar/Support/References/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa | \
gzip -c > ${outName}
