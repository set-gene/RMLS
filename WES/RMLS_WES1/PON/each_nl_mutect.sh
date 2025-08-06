#!/bin/bash

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/NL_1_bqsr.bam -max-mnp-distance 0 -O ./NL_1_pon.vcf.gz >> ./NL_1_pon.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/NL_2_bqsr.bam -max-mnp-distance 0 -O ./NL_2_pon.vcf.gz >> ./NL_2_pon.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/NL_3_bqsr.bam -max-mnp-distance 0 -O ./NL_3_pon.vcf.gz >> ./NL_3_pon.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/NL_4_bqsr.bam -max-mnp-distance 0 -O ./NL_4_pon.vcf.gz >> ./NL_4_pon.log 2>&1

gatk Mutect2 -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -I ../BQSR/NL_5_bqsr.bam -max-mnp-distance 0 -O ./NL_5_pon.vcf.gz >> ./NL_5_pon.log 2>&1
