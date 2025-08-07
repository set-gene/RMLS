#!/bin/bash

cat /HDD8T2/jeeh99/0118_RMLS/sort/sample.txt | while read id
do
 echo "start featureCounts for ${id}" 'date'
 featureCounts -T 6 -p -t exon -g gene_name -a /HDD2T/jeeh9/Reference/RNA_reference/gencode.v35.annotation.gtf -o /HDD8T2/jeeh99/0118_RMLS/DEG/gene_featurecounts/${id}.txt /HDD8T2/jeeh99/0118_RMLS/sort/${id}_sort.bam 2>&1|tee /HDD8T/jeeh99/0118_RMLS/DEG/gene_featurecounts/${id}_featurecounts.log
 echo "end featureCounts for ${id}" 'date'
done
