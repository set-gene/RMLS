#!/bin/bash

##FUSION1
cat /HDD8T2/jeeh99/0118_RMLS/FUSION_INSPECTOR/sample.txt | while read id
do
FusionInspector --fusions /HDD8T2/jeeh99/0118_RMLS/FUSION/${id}/${id}.fusionlist --genome_lib_dir /HDD2T/jeeh9/Reference/RNA_reference/CTAT_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir --left_fq /HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_1_P.fq.gz --right_fq /HDD8T/jeeh99/0114_RMLS_RNA/trimmed_fq/Paired/${id}_2_P.fq.gz --out_prefix ${id} --vis -O /HDD8T2/jeeh99/0118_RMLS/FUSION_INSPECTOR/${id} --no_remove_dups --CPU 10
done

##FUSION2
cat /HDD8T2/jeeh99/0118_RMLS/FUSION_INSPECTOR/sample2.txt | while read id
do
FusionInspector --fusions /HDD8T2/jeeh99/0118_RMLS/FUSION/${id}/${id}.fusionlist --genome_lib_dir /HDD2T/jeeh9/Reference/RNA_reference/CTAT_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir --left_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R1_P.fq.gz --right_fq /HDD2T/jeeh9/RMLS_RNA/trimmed_fq/Paired/${id}-R_R2_P.fq.gz --out_prefix ${id} --vis -O /HDD8T2/jeeh99/0118_RMLS/FUSION_INSPECTOR/${id} --no_remove_dups --CPU 10
done
