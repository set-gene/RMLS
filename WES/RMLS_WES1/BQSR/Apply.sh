#!/bin/bash

gatk ApplyBQSR -I ../Dedup/NL_1_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./NL_1_bqsr.table -O ./NL_1_bqsr.bam >> ./NL_1_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/NL_2_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./NL_2_bqsr.table -O ./NL_2_bqsr.bam >> ./NL_2_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/NL_3_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./NL_3_bqsr.table -O ./NL_3_bqsr.bam >> ./NL_3_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/NL_4_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./NL_4_bqsr.table -O ./NL_4_bqsr.bam >> ./NL_4_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/NL_5_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./NL_5_bqsr.table -O ./NL_5_bqsr.bam >> ./NL_5_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/MLS_1_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./MLS_1_bqsr.table -O ./MLS_1_bqsr.bam >> ./MLS_1_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/MLS_2_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./MLS_2_bqsr.table -O ./MLS_2_bqsr.bam >> ./MLS_2_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/MLS_3_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./MLS_3_bqsr.table -O ./MLS_3_bqsr.bam >> ./MLS_3_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/MLS_4_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./MLS_4_bqsr.table -O ./MLS_4_bqsr.bam >> ./MLS_4_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/MLS_5_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./MLS_5_bqsr.table -O ./MLS_5_bqsr.bam >> ./MLS_5_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/RMLS_1_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./RMLS_1_bqsr.table -O ./RMLS_1_bqsr.bam >> ./RMLS_1_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/RMLS_2_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./RMLS_2_bqsr.table -O ./RRMLS_2_bqsr.bam >> ./RMLS_2_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/RMLS_3_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./RMLS_3_bqsr.table -O ./RMLS_3_bqsr.bam >> ./RMLS_3_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/RMLS_4_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./RMLS_4_bqsr.table -O ./RMLS_4_bqsr.bam >> ./RMLS_4_apply.log 2>&1

gatk ApplyBQSR -I ../Dedup/RMLS_5_dedup.bam -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --bqsr-recal-file ./RMLS_5_bqsr.table -O ./RMLS_5_bqsr.bam >> ./RMLS_5_apply.log 2>&1
