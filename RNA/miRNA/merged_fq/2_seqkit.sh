#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/sample_raw.txt | while read id
do
	fq=/HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/${id}_merg.fq.gz
	echo "start seqkit for ${id}" 'date'
	seqkit seq -m 18 $fq | seqkit seq -M 40 > /HDD8T/jeeh99/0114_RMLS_RNA/miRNA/merged_fq/seqkit/${id}_filt.fastq.gz
	echo "end seqkit for ${id}" 'date'
done
