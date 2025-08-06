#!/bin/bash

sequenza-utils bam2seqz -n ./NL_1_sort.bam -t ./MLS_1_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/MLS_1.seqz.gz >> ../Sequenza_re/MLS_1_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_2_sort.bam -t ./MLS_2_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/MLS_2.seqz.gz >> ../Sequenza_re/MLS_2_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_3_sort.bam -t ./MLS_3_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/MLS_3.seqz.gz >> ../Sequenza_re/MLS_3_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_4_sort.bam -t ./MLS_4_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/MLS_4.seqz.gz >> ../Sequenza_re/MLS_4_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_5_sort.bam -t ./MLS_5_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/MLS_5.seqz.gz >> ../Sequenza_re/MLS_5_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_1_sort.bam -t ./RMLS_1_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/RMLS_1.seqz.gz >> ../Sequenza_re/RMLS_1_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_2_sort.bam -t ./RMLS_2_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/RMLS_2.seqz.gz >> ../Sequenza_re/RMLS_2_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_3_sort.bam -t ./RMLS_3_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/RMLS_3.seqz.gz >> ../Sequenza_re/RMLS_3_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_4_sort.bam -t ./RMLS_4_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/RMLS_4.seqz.gz >> ../Sequenza_re/RMLS_4_sequenza.log 2>&1

sequenza-utils bam2seqz -n ./NL_5_sort.bam -t ./RMLS_5_sort.bam --fasta /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -gc /HDD8T/jeeh99/GATK_REF/FASTA/hg38.gc50Base.wig.gz -o ../Sequenza_re/RMLS_5.seqz.gz >> ../Sequenza_re/RMLS_5_sequenza.log 2>&1
