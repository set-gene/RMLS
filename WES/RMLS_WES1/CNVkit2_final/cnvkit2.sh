#!/bin/bash

cnvkit.py batch -p 8 ../BQSR/MLS_1_bqsr.bam --normal ../BQSR/NL_1_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./MLS_1_reference.cnn --output-dir ./MLS_1 --drop-low-coverage --scatter --diagram >> ./MLS_1_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/MLS_2_bqsr.bam --normal ../BQSR/NL_2_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./MLS_2_reference.cnn --output-dir ./MLS_2 --drop-low-coverage --scatter --diagram >> ./MLS_2_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/MLS_3_bqsr.bam --normal ../BQSR/NL_3_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./MLS_3_reference.cnn --output-dir ./MLS_3 --drop-low-coverage --scatter --diagram >> ./MLS_3_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/MLS_4_bqsr.bam --normal ../BQSR/NL_4_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./MLS_4_reference.cnn --output-dir ./MLS_4 --drop-low-coverage --scatter --diagram >> ./MLS_4_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/MLS_5_bqsr.bam --normal ../BQSR/NL_5_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./MLS_5_reference.cnn --output-dir ./MLS_5 --drop-low-coverage --scatter --diagram >> ./MLS_5_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/RMLS_1_bqsr.bam --normal ../BQSR/NL_1_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./RMLS_1_reference.cnn --output-dir ./RMLS_1 --drop-low-coverage --scatter --diagram >> ./RMLS_1_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/RMLS_2_bqsr.bam --normal ../BQSR/NL_2_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./RMLS_2_reference.cnn --output-dir ./RMLS_2 --drop-low-coverage --scatter --diagram >> ./RMLS_2_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/RMLS_3_bqsr.bam --normal ../BQSR/NL_3_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./RMLS_3_reference.cnn --output-dir ./RMLS_3 --drop-low-coverage --scatter --diagram >> ./RMLS_3_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/RMLS_4_bqsr.bam --normal ../BQSR/NL_4_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./RMLS_4_reference.cnn --output-dir ./RMLS_4 --drop-low-coverage --scatter --diagram >> ./RMLS_4_cnvkit.log 2>&1

cnvkit.py batch -p 8 ../BQSR/RMLS_5_bqsr.bam --normal ../BQSR/NL_5_bqsr.bam --targets /HDD8T/jeeh99/GATK_REF/hg38.exon.bed --annotate /HDD8T/jeeh99/GATK_REF/GRCh38_gencode.v27.refFlat.txt --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta --method hybrid --output-reference ./RMLS_5_reference.cnn --output-dir ./RMLS_5 --drop-low-coverage --scatter --diagram >> ./RMLS_5_cnvkit.log 2>&1
