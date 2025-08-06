#!/bin/bash
cat ./sample.txt | while read id
do
echo "start ANNOVAR for ${id}" 'date'
~/tools/annovar/table_annovar.pl ./${id}.passed.somatic.vcf.gz ~/tools/annovar/humandb/ -buildver hg38 -out ../Annovar/${id} -remove -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp30a,cosmic96_coding,cosmic96_noncoding -operation gx,r,f,f,f,f,f -nastring . -vcfinput -thread 16 >> ../Annovar/${id}_annovar.log 2>&1
echo "end ANNOVAR for ${id}" 'date'
done
