#!/bin/bash

cat /HDD2T/jeeh9/RMLS_RNA/sort/sample.txt | while read id
do
 echo "start htseq-count for ${id}" 'date'
 htseq-count -f bam -r name -s no -a 10 -t gene -i gene_id -m union /HDD2T/jeeh9/RMLS_RNA/sort/${id}_sort.bam /HDD2T/jeeh9/Reference/RNA_reference/gencode.v35.annotation.gtf > /HDD2T/jeeh9/RMLS_RNA/DEG/htseq_gen/${id}.txt  
 echo "end htseq-count for ${id}" 'date'
done
