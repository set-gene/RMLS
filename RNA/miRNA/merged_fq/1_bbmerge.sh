#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/sample_raw.txt | while read id
do
	fq1=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_1_P.fq.gz
	fq2=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_2_P.fq.gz
	echo "start bbmerge for ${id}" 'date'
	bbmerge.sh in1=$fq1 in2=$fq2 out=/HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/${id}_merg.fq.gz 2>&1|tee /HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/${id}_merged.log
	echo "end bbmerge for ${id}" 'date'
done
