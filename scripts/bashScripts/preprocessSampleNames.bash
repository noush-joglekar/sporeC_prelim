#!/bin/bash
## Takes in a sample info text file and renames all the files in human readable format
## specified by input file

sampleInfo=$1;

while read line ; do \
sample=$(echo $line | awk '{sample="GM12878"_$1"_Rep"$2"_Read"$3".fastq.gz"; print sample}') ; \
id=$(echo $line | awk '{print $4}'); \
echo $sample; mv fastqs/$id.fastq.gz fastqs/$sample; done < $sampleInfo | tail -n+2
