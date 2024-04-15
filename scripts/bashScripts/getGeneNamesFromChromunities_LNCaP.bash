#!/bin/bash
# By Anoushka Joglekar
# 04.2023, modified 03.2024

## Takes in a fragments.csv file as outputted by Marcin's Pore-C pipeline
## and per read with cardinality >=2 it outputs gene names. These will then
## be fed into a script that looks at splicing efficiency for these gene sets.
## Can make it parallelizable to run per chromosome depending on run time.

## For multiple files I am using a batch ID in case of repeated read IDs

## New PoreC output is inconsistent with old one so making changes
## No way of getting read length and quality, which was "rl" in OG script

module load bedtools

fragFile=$1;
batch=$2;
outputDir=$3;

if [ ! -d "$outputDir" ]; then
mkdir $outputDir
fi

zcat $fragFile | awk -v batch=$batch '{FS=","; OFS="\t"; if($7=="True" && !($2~"_|EBV") ) {count[$1]++; \
if(coord[$1]){coord[$1]=coord[$1]"\n"$2"\t"$3"\t"$4"\tRead:"$1}else{coord[$1]=$2"\t"$3"\t"$4"\tRead:"$1}} }END \
{for(c in count) {if(count[c]>=2) {split(coord[c],entries,"\n"); \
for(en in entries) {print entries[en]"_Card:"count[c]"_Batch:"batch} }} }' > ${outputDir}/fragFile.bed

bedtools sort -i ${outputDir}/fragFile.bed > ${outputDir}/fragFile_sorted.bed
bedtools closest -d -a ${outputDir}/fragFile_sorted.bed -b geneCoordsAndNames_sorted.bed | tr -d "\"|;" > ${outputDir}/fragFile_sorted_wClosestGene.bed
cat ${outputDir}/fragFile_sorted_wClosestGene.bed | awk '!seen[$1"_"$2"_"$3"_"$4]++' > ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed
gzip ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed

time zcat ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed.gz | awk '{split($4,r,/:|_/); OFS="\t"; \
print r[5]":"r[6]"_"r[1]":"r[2],r[4],$9,$10,$11,$11}' | \
gzip -c > ${outputDir}/readChromunities_wClosestGene_Batch${batch}.gz

rm ${outputDir}/fragFile.bed ${outputDir}/fragFile_sorted.bed ${outputDir}/fragFile_sorted_wClosestGene.bed