#!/bin/bash

gatk CreateSomaticPanelOfNormals -R /HDD8T/jeeh99/GATK_REF/FASTA/Homo_sapiens_assembly38.fasta -V gendb://pon_db -O ./NL_all_pon.vcf.gz >> ./mk_pon.log 2>&1
