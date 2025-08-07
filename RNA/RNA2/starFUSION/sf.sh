#!/bin/bash
##FUSION1
cat /HDD8T2/jeeh99/0118_RMLS/FUSION/sample1.txt | while read id
do
STAR-Fusion --genome_lib_dir /HDD2T/jeeh9/Reference/RNA_reference/CTAT_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir --left_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R1_P.fq.gz --right_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R2_P.fq.gz --output_dir /HDD8T2/jeeh99/0118_RMLS/FUSION/${id} --no_remove_dups --CPU 10
done
##FUSION2
cat /HDD8T2/jeeh99/0118_RMLS/FUSION/sample2.txt | while read id
do
STAR-Fusion --genome_lib_dir /HDD2T/jeeh9/Reference/RNA_reference/CTAT_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir --left_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R1_P.fq.gz --right_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R2_P.fq.gz --output_dir /HDD8T2/jeeh99/0118_RMLS/FUSION/${id} --no_remove_dups --CPU 10
done
