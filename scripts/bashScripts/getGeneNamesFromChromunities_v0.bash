#!/bin/bash
# By Anoushka Joglekar
# 04.2023

## Takes in a fragments.csv file as outputted by Marcin's Pore-C pipeline
## and per read with cardinality >=2 it outputs gene names. These will then
## be fed into a script that looks at splicing efficiency for these gene sets.
## Can make it parallelizable to run per chromosome depending on run time.

module load bedtools

fragFile=$1;

zcat $fragFile | awk '{FS=","; OFS="\t"; if($18=="True" && !($4~"_|EBV") ) {count[$1]++; \
if(coord[$1]){coord[$1]=coord[$1]"\n"$4"\t"$5"\t"$6"\tRead:"$1}else{coord[$1]=$4"\t"$5"\t"$6"\tRead:"$1}} }END \
{for(c in count) {if(count[c]>=2) {split(coord[c],entries,"\n"); for(en in entries) {print entries[en]"_Card:"count[c]} }} }' > fragFile.bed

bedtools sort -i fragFile.bed > fragFile_sorted.bed
bedtools closest -d -a fragFile_sorted.bed -b geneCoordsAndNames_sorted.bed | tr -d "\"|;" > fragFile_sorted_wClosestGene.bed
cat fragFile_sorted_wClosestGene.bed | awk '!seen[$1"_"$2"_"$3"_"$4]++' > x ; mv x fragFile_sorted_wClosestGene.bed
gzip fragFile_sorted_wClosestGene.bed

time zcat fragFile_sorted_wClosestGene.bed.gz | awk '{split($4,r,/:|_/); OFS="\t"; print r[1]":"r[2],r[4],$9,$10,$11,$12}' | \
gzip -c > readChromunities_wClosestGene.gz
