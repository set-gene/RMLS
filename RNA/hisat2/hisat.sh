#!/bin/bash

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/NL_sample.txt | while read id
do
        fq1=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_1_P.fq.gz
        fq2=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_2_P.fq.gz
        echo "start hisat2 for ${id}" 'date'
        hisat2 -p 8 -x /HDD2T/jeeh9/Reference/RNA_reference/gencode_genome/gencode_index -1 ${fq1} -2 ${fq2} -S /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.sam >> /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.log 2>&1
        echo "end hisat2 for ${id}" 'date'
done

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/MLS_sample.txt | while read id
do
	fq1=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_1_P.fq.gz
	fq2=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_2_P.fq.gz
	echo "start hisat2 for ${id}" 'date'
	hisat2 -p 8 -x /HDD2T/jeeh9/Reference/RNA_reference/gencode_genome/gencode_index -1 ${fq1} -2 ${fq2} -S /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.sam >> /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.log 2>&1
	echo "end hisat2 for ${id}" 'date'
done

cat /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/RMLS_sample.txt | while read id
do
        fq1=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_1_P.fq.gz
        fq2=/HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_2_P.fq.gz
        echo "start hisat2 for ${id}" 'date'
        hisat2 -p 8 -x /HDD2T/jeeh9/Reference/RNA_reference/gencode_genome/gencode_index -1 ${fq1} -2 ${fq2} -S /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.sam >> /HDD8T/jeeh99/0114_RMLS_RNA/hisat2/${id}.log 2>&1
        echo "end hisat2 for ${id}" 'date'
done
