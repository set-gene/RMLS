#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/NL_sample.txt | while read id 
do
echo "start samtools sort name for ${id}" 'date'
samtools sort -@ 10 /HDD8T/jeeh99/0114_RMLS_RNA/BAM/${id}.bam -o /HDD8T/jeeh99/0114_RMLS_RNA/sort/${id}_sort.bam
echo "end samtools sort name for ${id}" 'date'
done

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/MLS_sample.txt | while read id 
do
echo "start samtools sort name for ${id}" 'date'
samtools sort -@ 10 /HDD8T/jeeh99/0114_RMLS_RNA/BAM/${id}.bam -o /HDD8T/jeeh99/0114_RMLS_RNA/sort/${id}_sort.bam
echo "end samtools sort name for ${id}" 'date'
done

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/RMLS_sample.txt | while read id 
do
echo "start samtools sort name for ${id}" 'date'
samtools sort -@ 10 /HDD8T/jeeh99/0114_RMLS_RNA/BAM/${id}.bam -o /HDD8T/jeeh99/0114_RMLS_RNA/sort/${id}_sort.bam
echo "end samtools sort name for ${id}" 'date'
done
