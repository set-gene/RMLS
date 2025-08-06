#!/bin/bash

gatk BaseRecalibrator -I ../Dedup/NL_1_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./NL_1_bqsr.table >> ./NL_1.log 2>&1
	
gatk BaseRecalibrator -I ../Dedup/NL_2_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./NL_2_bqsr.table >> ./NL_2.log 2>&1

gatk BaseRecalibrator -I ../Dedup/NL_3_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./NL_3_bqsr.table >> ./NL_3.log 2>&1

gatk BaseRecalibrator -I ../Dedup/NL_4_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./NL_4_bqsr.table >> ./NL_4.log 2>&1

gatk BaseRecalibrator -I ../Dedup/NL_5_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./NL_5_bqsr.table >> ./NL_5.log 2>&1

gatk BaseRecalibrator -I ../Dedup/MLS_1_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./MLS_1_bqsr.table >> ./MLS_1.log 2>&1
	
gatk BaseRecalibrator -I ../Dedup/MLS_2_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./MLS_2_bqsr.table >> ./MLS_2.log 2>&1

gatk BaseRecalibrator -I ../Dedup/MLS_3_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./MLS_3_bqsr.table >> ./MLS_3.log 2>&1

gatk BaseRecalibrator -I ../Dedup/MLS_4_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./MLS_4_bqsr.table >> ./MLS_4.log 2>&1

gatk BaseRecalibrator -I ../Dedup/MLS_5_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./ MLS_5_bqsr.table >> ./MLS_5.log 2>&1

gatk BaseRecalibrator -I ../Dedup/RMLS_1_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./RMLS_1_bqsr.table >> ./RMLS_1.log 2>&1
	
gatk BaseRecalibrator -I ../Dedup/RMLS_2_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./RMLS_2_bqsr.table >> ./RMLS_2.log 2>&1

gatk BaseRecalibrator -I ../Dedup/RMLS_3_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./RMLS_3_bqsr.table >> ./RMLS_3.log 2>&1

gatk BaseRecalibrator -I ../Dedup/RMLS_4_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./RMLS_4_bqsr.table >> ./RMLS_4.log 2>&1

gatk BaseRecalibrator -I ../Dedup/RMLS_5_dedup.bam --known-sites /HDD8T/jeeh99/GATK_REF/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /HDD8T/jeeh99/GATK_REF/1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /HDD8T/jeeh99/GATK_REF/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -O ./RMLS_5_bqsr.table >> ./RMLS_5.log 2>&1
