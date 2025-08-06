#!/bin/bash
cat /HDD8T3/jeeh99/RMLS_WES1/sort/RMLS.txt | while read id

do 

samtools index -b -@ 10 ${id}_sort.bam ${id}_sort.bam.bai

samtools flagstat ${id}_sort.bam > ${id}.flagstat

samtools idxstats ${id}_sort.bam > ${id}_sort.bam.idx

done  
