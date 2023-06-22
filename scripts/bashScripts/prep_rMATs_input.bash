#!/bin/bash

baseName=$1;
outName=$2;

#mkdir ${baseName};
a=$(ls -m ${baseName}*/STAR_alignment/${baseName}*.Aligned.out.WithReadGroup.sorted.markdup.bam); \
echo $a | tr -d ' ' > $outName
