#!/bin/bash

sequenza-utils seqz_binning --seqz MLS_1.seqz.gz -w 50 -o MLS_1.small.seqz.gz >> ./MLS_1_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz MLS_2.seqz.gz -w 50 -o MLS_2.small.seqz.gz >> ./MLS_2_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz MLS_3.seqz.gz -w 50 -o MLS_3.small.seqz.gz >> ./MLS_3_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz MLS_4.seqz.gz -w 50 -o MLS_4.small.seqz.gz >> ./MLS_4_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz MLS_5.seqz.gz -w 50 -o MLS_5.small.seqz.gz >> ./MLS_5_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz RMLS_1.seqz.gz -w 50 -o RMLS_1.small.seqz.gz >> ./RMLS_1_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz RMLS_2.seqz.gz -w 50 -o RMLS_2.small.seqz.gz >> ./RMLS_2_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz RMLS_3.seqz.gz -w 50 -o RMLS_3.small.seqz.gz >> ./RMLS_3_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz RMLS_4.seqz.gz -w 50 -o RMLS_4.small.seqz.gz >> ./RMLS_4_seqbin.log 2>&1

sequenza-utils seqz_binning --seqz RMLS_5.seqz.gz -w 50 -o RMLS_5.small.seqz.gz >> ./RMLS_5_seqbin.log 2>&1

