#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/NL_sample.txt | while read id
do
echo "start samtools view for ${id}" 'date'
samtools view -Sb /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.sam > /HDD8T/jeeh99/0114_RMLS_RNA/BAM/${id}.bam
echo "end samtools view for ${id}" 'date'
done
cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/MLS_sample.txt | while read id
do
echo "start samtools view for ${id}" 'date'
samtools view -Sb /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.sam > /HDD8T/jeeh99/0114_RMLS_RNA/BAM/${id}.bam
echo "end samtools view for ${id}" 'date'
done

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/RMLS_sample.txt | while read id
do
echo "start samtools view for ${id}" 'date'
samtools view -Sb /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.sam > /HDD8T/jeeh99/0114_RMLS_RNA/BAM/${id}.bam
echo "end samtools view for ${id}" 'date'
done
