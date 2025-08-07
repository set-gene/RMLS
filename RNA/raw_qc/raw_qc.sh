#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/raw_fq/sample_raw.txt | while read id
do
	fastqc --outdir ./raw_qc/ --threads 16 ./${id}*.fastq.gz >> ./raw_qc/${id}_fastqc.log 2>&1 
done 
