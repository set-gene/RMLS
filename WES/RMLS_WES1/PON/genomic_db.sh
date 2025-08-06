#!/bin/bash

gatk GenomicsDBImport -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -L chr1 -L chr2 -L chr3  -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY --genomicsdb-workspace-path pon_db -V ./NL_1_pon.vcf.gz -V ./NL_2_pon.vcf.gz -V ./NL_3_pon.vcf.gz -V ./NL_4_pon.vcf.gz -V ./NL_5_pon.vcf.gz -V ./NL_6_pon.vcf.gz -V ./NL_7_pon.vcf.gz -V ./NL_8_pon.vcf.gz -V ./NL_9_pon.vcf.gz >> ./db_pon.log 2>&1
