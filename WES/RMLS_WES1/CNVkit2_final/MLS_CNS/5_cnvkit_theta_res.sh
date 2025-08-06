#!/bin/sh
cat config | while read line 
do 
  tumor=$(echo $line | cut -d ' ' -f 2)
  
    RunTHetA.py ./${tumor}_call.filt.interval_count \
    --TUMOR_FILE ./${tumor}_call.filt.tumor.snp_formatted.txt \
    --NORMAL_FILE ./${tumor}_call.filt.normal.snp_formatted.txt \
    --BAF --NUM_PROCESSES `nproc` --FORCE
done
