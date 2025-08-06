#!/bin/bash

bwa mem -t 10 -M -R "@RG\tID:MLS_1\tSM:MLS_1\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_1_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_1_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_1.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_1.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:MLS_2\tSM:MLS_2\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_2_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_2_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_2.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_2.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:MLS_3\tSM:MLS_3\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_3_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_3_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_3.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_3.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:MLS_4\tSM:MLS_4\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_4_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_4_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_4.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_4.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:MLS_5\tSM:MLS_5\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_5_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/MLS_5_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_5.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/MLS_5.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:NL_1\tSM:NL_1\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_1_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_1_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_1.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_1.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:NL_2\tSM:NL_2\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_2_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_2_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_2.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_2.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:NL_3\tSM:NL_3\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_3_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_3_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_3.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_3.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:NL_4\tSM:NL_4\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_4_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_4_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_4.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_4.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:NL_5\tSM:NL_5\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_5_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/NL_5_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_5.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/NL_5.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:RMLS_1\tSM:RMLS_1\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_1_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_1_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_1.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_1.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:RMLS_2\tSM:RMLS_2\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_2_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_2_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_2.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_2.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:RMLS_3\tSM:RMLS_3\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_3_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_3_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_3.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_3.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:RMLS_4\tSM:RMLS_4\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_4_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_4_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_4.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_4.log 2>&1

bwa mem -t 10 -M -R "@RG\tID:RMLS_5\tSM:RMLS_5\tLB:WES\tPL:ILLUMINA" /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_5_1_P.fq.gz /HDD8T/jeeh99/RMLS_WES1/trimmed_fq/Paired/RMLS_5_2_P.fq.gz -o /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_5.sam >> /HDD8T3/jeeh99/RMLS_WES1/bwa/RMLS_5.log 2>&1
