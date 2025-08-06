 #!/usr/bin/env Rscript

library(argparse)
library(copynumber)
library(sequenza)

parser = ArgumentParser()
parser$add_argument('-i', '--input', required=TRUE) # mysample.small.seqz.gz

# result/mysample
args = parser$parse_args()

chromosome_list<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

extracted = sequenza.extract(args$input,assembly = "hg38",chromosome.list=chromosome_list)

CP = sequenza.fit(extracted)

confint <- get.ci(CP)

purity <- confint$max.cellularity

print(purity)
