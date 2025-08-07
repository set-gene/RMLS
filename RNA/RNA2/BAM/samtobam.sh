#!/bin/bash

cat /HDD8T2/jeeh99/0118_RMLS/hisat2/mls.txt | while read id
do
echo "start samtools view for ${id}" 'date'
samtools view -Sb /HDD8T2/jeeh99/0118_RMLS/hisat2/${id}.sam > /HDD8T2/jeeh99/0118_RMLS/BAM/${id}.bam
echo "end samtools view for ${id}" 'date'
done
