#!/bin/bash

cat config_all | while read line 
do 
  tumor=$(echo $line | cut -d ' ' -f 2)

cnvkit.py export seg ./${tumor}_call.filt.cns\
  -o ./${tumor}.filt.seg   
done
