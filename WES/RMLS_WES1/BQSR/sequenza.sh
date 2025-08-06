#!/bin/bash

sequenza-utils bam2seqz -n ./NL_1_bqsr.bam -t ./MLS_1_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/MLS_1.seqz.gz >> ../Sequenza/MLS_1_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_2_bqsr.bam -t ./MLS_2_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/MLS_2.seqz.gz >> ../Sequenza/MLS_2_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_3_bqsr.bam -t ./MLS_3_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/MLS_3.seqz.gz >> ../Sequenza/MLS_3_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_4_bqsr.bam -t ./MLS_4_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/MLS_4.seqz.gz >> ../Sequenza/MLS_4_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_5_bqsr.bam -t ./MLS_5_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/MLS_5.seqz.gz >> ../Sequenza/MLS_5_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_1_bqsr.bam -t ./RMLS_1_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/RMLS_1.seqz.gz >> ../Sequenza/RMLS_1_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_2_bqsr.bam -t ./RMLS_2_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/RMLS_2.seqz.gz >> ../Sequenza/RMLS_2_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_3_bqsr.bam -t ./RMLS_3_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/RMLS_3.seqz.gz >> ../Sequenza/RMLS_3_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_4_bqsr.bam -t ./RMLS_4_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/RMLS_4.seqz.gz >> ../Sequenza/RMLS_4_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_5_bqsr.bam -t ./RMLS_5_bqsr.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza/RMLS_5.seqz.gz >> ../Sequenza/RMLS_5_sequenza.log 2>&1
