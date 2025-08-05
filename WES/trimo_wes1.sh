#!/bin/bash
cat /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/sample_raw.txt | while read id 
do
	fq1=/HDD8T/jeeh99/RMLS_WES1/raw_fq/file_${id}/${id}_1.fastq.gz
	fq2=/HDD8T/jeeh99/RMLS_WES1/raw_fq/file_${id}/${id}_2.fastq.gz
	echo "start trimmmomatic for ${id}" 'date'
	java -jar ~/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 10 ${fq1} ${fq2} /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/file_${id}_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Unpair/file_${id}_1_U.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/file_${id}_2_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Unpair/file_${id}_2_U.fq.gz ILLUMINACLIP:/home/jeeh99/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50 >> /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/file_${id}.log 2>&1
	echo "end trimmomatic for ${id}" 'date'
done
