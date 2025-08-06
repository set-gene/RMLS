#!/bin/bash
gunzip NLsamplepon_.vcf.gz 

cd /home/tools

java -jar picard.jar SortVcf I=/HDD2T/jeeh9/pon/NLPON/NLsamplepon_.vcf O=/HDD2T/jeeh9/pon/NLPON/NLsamplepon_sort.vcf
