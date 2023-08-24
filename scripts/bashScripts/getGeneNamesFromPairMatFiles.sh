#!/bin/bash
# By Anoushka Joglekar
# 04.2023

## Takes in a fragments.csv file as outputted by Marcin's Pore-C pipeline
## and per read with cardinality >=2 it outputs gene names. These will then
## be fed into a script that looks at splicing efficiency for these gene sets.
## Can make it parallelizable to run per chromosome depending on run time.

## For multiple files I am using a batch ID in case of repeated read IDs
## In python can assign new IDs as well.
## also an output dir needs to be specified

module load bedtools

fragFile=$1;
batch=$2;
outputDir=$3;

if [ ! -d "$outputDir" ]; then
mkdir $outputDir
fi

zcat $fragFile | awk -v batch=$batch '{FS=","; OFS="\t"; {count[$1]++; rl[$1]=$9"\t"$12; \
if(coord[$1]){coord[$1]=coord[$1]"\n"$4"\t"$5"\t"$6"\tRead:"$1}else{coord[$1]=$4"\t"$5"\t"$6"\tRead:"$1}} }END \
{for(c in count) {if(count[c]>=2) {split(coord[c],entries,"\n"); \
for(en in entries) {print entries[en]"_Card:"count[c]"_Batch:"batch"\t"rl[c]} }} }' > fragFile.bed

bedtools sort -i fragFile.bed > fragFile_sorted.bed
bedtools closest -d -a fragFile_sorted.bed -b geneCoordsAndNames_sorted.bed | tr -d "\"|;" > fragFile_sorted_wClosestGene.bed
cat fragFile_sorted_wClosestGene.bed | awk '!seen[$1"_"$2"_"$3"_"$4]++' > ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed
gzip ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed

time zcat ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed.gz | awk '{split($4,r,/:|_/); OFS="\t"; \
print r[5]":"r[6]"_"r[1]":"r[2],r[4],$5,$6,$11,$12,$13,$14}' | \
gzip -c > ${outputDir}/readChromunities_wClosestGene_Batch${batch}.gz

rm fragFile.bed fragFile_sorted.bed fragFile_sorted_wClosestGene.bed
