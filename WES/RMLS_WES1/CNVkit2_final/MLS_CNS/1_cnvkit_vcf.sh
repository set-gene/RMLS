#!/bin/sh
cat config | while read line 
do 
  normal=$(echo $line | cut -d ' ' -f 1)
  tumor=$(echo $line | cut -d ' ' -f 2)

  cat ../../filt/${tumor}.somatic.filt.vcf | \
  awk 'BEGIN{FS="\t";OFS="\t"} {if($0~"^#") {print $0} else {if($7!="PASS"){$7="REJECT"} print $0}}'  > ./${tumor}_cnvkit.vcf
  sed -i '2 i ##FILTER=<ID=REJECT,Description="Not somatic due to normal call frequency or phred likelihoods: tumor: 35, normal 35.">' ./${tumor}_cnvkit.vcf
  sed -i "2 i ##PEDIGREE=<Derived=${tumor},Original=${normal}>" ./${tumor}_cnvkit.vcf
  
done
