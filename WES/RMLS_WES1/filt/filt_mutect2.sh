#!/bin/bash


cat ./sample.txt | while read id
do
gatk FilterMutectCalls -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -V ../Mutect2/${id}_somatic.vcf.gz -O ./${id}.somatic.filt.vcf
done

gunzip *.gz
#####################################pass filter
cat ./sample.txt | while read id
do
cat ./${id}.somatic.filt.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ${id}_pass.somatic.filt.vcf 
done

######################################select
cat ./sample.txt | while read id
do
gatk SelectVariants -select "DP>30" -V ./${id}_pass.somatic.filt.vcf -O ./${id}.passed.somatic.vcf.gz
done
