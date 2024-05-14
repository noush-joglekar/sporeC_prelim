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

## THINGS SEEM TO NOT MATCH IN QC. First column is random ID, does not correspond to read ID

module load bedtools

fragFile=$1;
batch=$2;
outputDir=$3;

if [ ! -d "$outputDir" ]; then
mkdir $outputDir
fi

zcat $fragFile | awk -v batch=$batch '{FS=","; OFS="\t"; if($6!=current) {id++; current=$6;} \
if($7=="True" && !($2~"_|EBV") ) {count[id]++; \
if(coord[id]){coord[id]=coord[id]"\n"$2"\t"$3"\t"$4"\tRead:"id}else{coord[id]=$2"\t"$3"\t"$4"\tRead:"id}} }END \
{for(c in count) {if(count[c]>=2) {split(coord[c],entries,"\n"); \
for(en in entries) {print entries[en]"_Card:"count[c]"_Batch:"batch} }} }' > ${outputDir}/fragFile_${batch}.bed

bedtools sort -i ${outputDir}/fragFile_${batch}.bed > ${outputDir}/fragFile_${batch}_sorted.bed
bedtools closest -d -a ${outputDir}/fragFile_${batch}_sorted.bed -b geneCoordsAndNames_sorted.bed | tr -d "\"|;" > ${outputDir}//fragFile_sorted_wClosestGene_tmp_${batch}.bed
cat ${outputDir}//fragFile_sorted_wClosestGene_tmp_${batch}.bed | awk '!seen[$1"_"$2"_"$3"_"$4]++' > ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed
gzip ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed

time zcat ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed.gz | awk '{split($4,r,/:|_/); OFS="\t"; \
print r[5]":"r[6]"_"r[1]":"r[2],r[4],$9,$10,$11,$11}' | \
gzip -c > ${outputDir}/readChromunities_wClosestGene_Batch${batch}.gz

rm ${outputDir}/fragFile_${batch}.bed ${outputDir}/fragFile_${batch}_sorted.bed ${outputDir}//fragFile_sorted_wClosestGene_tmp_${batch}.bed