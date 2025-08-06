#!/bin/sh
cat config | while read line 
do 
  tumor=$(echo $line | cut -d ' ' -f 2)
  
cnvkit.py import-theta ./${tumor}_call.filt.cns ./${tumor}_call.BEST.results
    
done
