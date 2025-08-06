#!/bin/bash

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/NL_1_sort.bam O=./NL_1_dedup.bam M=./NL_1_dedup.metrics.txt

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/NL_2_sort.bam O=./NL_2_dedup.bam M=./NL_2_dedup.metrics.txt

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/NL_3_sort.bam O=./NL_3_dedup.bam M=./NL_3_dedup.metrics.txt

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/NL_4_sort.bam O=./NL_4_dedup.bam M=./NL_4_dedup.metrics.txt

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/NL_5_sort.bam O=./NL_5_dedup.bam M=./NL_5_dedup.metrics.txt

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/MLS_1_sort.bam O=./MLS_1_dedup.bam M=./MLS_1_dedup.metrics.txt >> ./MLS_1.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/MLS_2_sort.bam O=./MLS_2_dedup.bam M=./MLS_2_dedup.metrics.txt >> ./MLS_2.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/MLS_3_sort.bam O=./MLS_3_dedup.bam M=./MLS_3_dedup.metrics.txt >> ./MLS_3.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/MLS_4_sort.bam O=./MLS_4_dedup.bam M=./MLS_4_dedup.metrics.txt >> ./MLS_4.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/MLS_5_sort.bam O=./MLS_5_dedup.bam M=./MLS_5_dedup.metrics.txt >> ./MLS_5.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/RMLS_1_sort.bam O=./RMLS_1_dedup.bam M=./RMLS_1_dedup.metrics.txt >> ./RMLS_1.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/RMLS_2_sort.bam O=./RMLS_2_dedup.bam M=./RMLS_2_dedup.metrics.txt >> ./RMLS_2.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/RMLS_3_sort.bam O=./RMLS_3_dedup.bam M=./RMLS_3_dedup.metrics.txt >> ./RMLS_3.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/RMLS_4_sort.bam O=./RMLS_4_dedup.bam M=./ RMLS_4_dedup.metrics.txt >> ./RMLS_4.log 2>&1

java -jar ~/tools/picard.jar MarkDuplicates I=../sort/RMLS_5_sort.bam O=./RMLS_5_dedup.bam M=./RMLS_5_dedup.metrics.txt >> ./RMLS_5.log 2>&1
