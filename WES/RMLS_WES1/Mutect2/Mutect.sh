#!/bin/bash

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/MLS_1_bqsr.bam -I ../BQSR/NL_1_bqsr.bam -normal NL_1 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./MLS_1_somatic.vcf.gz >> ./mls_1_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/MLS_2_bqsr.bam -I ../BQSR/NL_2_bqsr.bam -normal NL_2 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./MLS_2_somatic.vcf.gz >> ./mls_2_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/MLS_3_bqsr.bam -I ../BQSR/NL_3_bqsr.bam -normal NL_3 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./MLS_3_somatic.vcf.gz >> ./mls_3_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/MLS_4_bqsr.bam -I ../BQSR/NL_4_bqsr.bam -normal NL_4 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./MLS_4_somatic.vcf.gz >> ./mls_4_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/MLS_5_bqsr.bam -I ../BQSR/NL_5_bqsr.bam -normal NL_5 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./MLS_5_somatic.vcf.gz >> ./mls_5_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/RMLS_1_bqsr.bam -I ../BQSR/NL_1_bqsr.bam -normal NL_1 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./RMLS_1_somatic.vcf.gz >> ./rmls_1_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/RMLS_2_bqsr.bam -I ../BQSR/NL_2_bqsr.bam -normal NL_2 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./RMLS_2_somatic.vcf.gz >> ./rmls_2_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/RMLS_3_bqsr.bam -I ../BQSR/NL_3_bqsr.bam -normal NL_3 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./RMLS_3_somatic.vcf.gz >> ./rmls_3_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/RMLS_4_bqsr.bam -I ../BQSR/NL_4_bqsr.bam -normal NL_4 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./RMLS_4_somatic.vcf.gz >> ./rmls_4_vcf.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/RMLS_5_bqsr.bam -I ../BQSR/NL_5_bqsr.bam -normal NL_5 --germline-resource /HDD8T/jeeh99/GATK_REF/af-only-gnomad.hg38.vcf.gz --panel-of-normals ../PON_RE/NL_all_pon_sort.vcf -O ./RMLS_5_somatic.vcf.gz >> ./rmls_5_vcf.log 2>&1
