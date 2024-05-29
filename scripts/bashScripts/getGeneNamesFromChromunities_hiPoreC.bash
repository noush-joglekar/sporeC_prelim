#!/bin/bash
# By Anoushka Joglekar
# 05.2024

## Yet another format -__-

module load bedtools

fragFile=$1;
batch=$2;
outputDir=$3;

if [ ! -d "$outputDir" ]; then
mkdir $outputDir
fi

zcat $fragFile | awk -v batch=$batch '{FS=","; if($1!=current) {id++; current=$1;} \
if(!($2~"_|EBV") ) {count[id]++; \
if(coord[id]){coord[id]=coord[id]"\n"$6"\t"$8"\t"$9"\tRead:"id}else{coord[id]=$6"\t"$8"\t"$9"\tRead:"id}} }END \
{for(c in count) {if(count[c]>=2) {split(coord[c],entries,"\n"); \
for(en in entries) {print entries[en]"_Card:"count[c]"_Batch:"batch} }} }' > ${outputDir}/fragFile_${batch}.bed

bedtools sort -i ${outputDir}/fragFile_${batch}.bed > ${outputDir}/fragFile_${batch}_sorted.bed
bedtools closest -d -a ${outputDir}/fragFile_${batch}_sorted.bed -b geneCoordsAndNames_sorted.bed | tr -d "\"|;" > ${outputDir}/fragFile_sorted_wClosestGene_tmp_${batch}.bed

cat ${outputDir}/fragFile_sorted_wClosestGene_tmp_${batch}.bed | awk '!seen[$1"_"$2"_"$3"_"$4]++' > ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed
gzip ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed

time zcat ${outputDir}/fragFile_sorted_wClosestGene_Batch${batch}.bed.gz | awk '{split($4,r,/:|_/); OFS="\t"; \
print r[5]":"r[6]"_"r[1]":"r[2],r[4],$9,$10,$11,$11}' | \
gzip -c > ${outputDir}/readChromunities_wClosestGene_Batch${batch}.gz

#rm ${outputDir}/fragFile_${batch}.bed ${outputDir}/fragFile_${batch}_sorted.bed ${outputDir}/fragFile_sorted_wClosestGene_tmp_${batch}.bed
