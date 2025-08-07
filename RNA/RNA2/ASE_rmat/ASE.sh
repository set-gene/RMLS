# analysis of Alternative splicing events

rmats.py --b1 ./ASE_M_R/MLS.txt --b2 ./ASE_M_R/RMLS.txt --gtf /HDD2T/jeeh9/Reference/RNA_reference/gencode.v35.annotation.gtf  --od ./ASE_M_R/output -t paired --readLength 101 --cstat 0.0001 --nthread 10 

rmats.py --b1 ./ASE_N_M/NL.txt --b2 ./ASE_N_M/MLS.txt --gtf /HDD2T/jeeh
9/Reference/RNA_reference/gencode.v35.annotation.gtf  --od ./ASE_N_M/output -t p
aired --readLength 101 --cstat 0.0001 --nthread 10

rmats.py --b1 ./ASE_N_R/NL.txt --b2 ./ASE_N_R/RMLS.txt --gtf /HDD2T/jeeh
9/Reference/RNA_reference/gencode.v35.annotation.gtf  --od ./ASE_N_R/output -t p
aired --readLength 101 --cstat 0.0001 --nthread 10
