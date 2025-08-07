#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/sample_raw.txt | while read id
do
	fq=/HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/fastqq/${id}_filt.fastq
	echo "start fastx for ${id}" 'date'
	fastx_trimmer -v -f 1 -l 27 -i $fq -Q 33 -z -o /HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/fast_clean/${id}_clean.fq.gz
	echo "end fastx for ${id}" 'date'
done
