#!/bin/bash

cat /HDD8T2/jeeh99/0118_RMLS/hisat2/sample.txt | while read id 
do
echo "start samtools sort name for ${id}" 'date'
samtools sort -@ 10 /HDD8T2/jeeh99/0118_RMLS/BAM/${id}.bam -o /HDD8T2/jeeh99/0118_RMLS/sort/${id}_sort.bam
echo "end samtools sort name for ${id}" 'date'
done
