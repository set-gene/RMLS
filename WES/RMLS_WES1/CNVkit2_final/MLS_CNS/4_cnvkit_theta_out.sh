#!/bin/sh
cat config_all | while read line 
do 
  tumor=$(echo $line | cut -d ' ' -f 2)
  
  cnvkit.py export theta ./${tumor}_call.filt.cns -r ../${tumor}_reference.cnn -v ./${tumor}_cnvkit.vcf
  
done
