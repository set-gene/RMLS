#!/bin/bash
cat /HDD2T/jeeh9/RMLS_RNA/sort/sample.txt | while read id
do
 echo "start RSEM for ${id}" 'date'
rsem-calculate-expression --paired-end --bam --estimate-rspd --append-names --no-bam-output -p 6 /HDD2T/jeeh9/RMLS_RNA/sort/${id}_sort.bam /HDD2T/jeeh9/Reference/RNA_Reference/ensembl/RSEM_idx /HDD2T/jeeh9/RMLS_RNA/isoform/RSEM/${id}  >> /HDD2T/jeeh9/RMLS_RNA/isoform/RSEM/${id}_RSEM.log 2>&1  
 echo "end RSEM for ${id}" 'date'
done
