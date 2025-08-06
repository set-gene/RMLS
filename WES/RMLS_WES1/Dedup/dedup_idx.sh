cat /HDD8T3/jeeh99/RMLS_WES1/sort/NL.txt | while read id

do 
samtools index -b -@ 10 ./${id}_dedup.bam ./${id}_dedup.bam.bai
done  
