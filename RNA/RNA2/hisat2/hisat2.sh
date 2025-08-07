#!/bin/bash

cat /HDD8T2/jeeh99/0118_RMLS/hisat2/mls.txt | while read id
do
        fq1=/HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}_R1_P.fq.gz
        fq2=/HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}_R2_P.fq.gz
        echo "start hisat2 for ${id}" 'date'
        hisat2 -p 8 -x /HDD2T/jeeh9/Reference/RNA_reference/gencode_genome/gencode_index -1 ${fq1} -2 ${fq2} -S /HDD8T2/jeeh99/0118_RMLS/hisat2/${id}.sam >> /HDD8T2/jeeh99/0118_RMLS/hisat2/${id}.log 2>&1
        echo "end hisat2 for ${id}" 'date'
done
