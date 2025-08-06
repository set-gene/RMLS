#!/bin/bash

cat /HDD8T/jeeh99/RMLS_WES1/raw_fq/sample_raw.txt | while read id
do
	fastqc --outdir ./raw_qc/ --threads 16 ./file_${id}/${id}*.fastq.gz >> ./raw_qc/${id}_fastqc.log 2>&1 
done 
