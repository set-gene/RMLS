MLS <- read.table("RMLS_9.hg38_multianno.txt",sep = "\t",header = TRUE)


MLS2  <-MLS[,c(7,60:70)]

colnames(MLS2 ) <- c("Gene","CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NL","RMLS")

sample_id <- "RMLS_9"


somatic <- data.frame(mutation_id=paste(MLS2$CHROM,MLS2$POS,MLS2$Gene, sep=":"), 
                      sample_id = sample_id,
                      chr=MLS2$CHROM,
                      start=MLS2$POS,
                      end=MLS2$POS+nchar(MLS2$REF),
                      ref_counts=sapply(
                        sapply(sapply(strsplit(MLS2$RMLS,":"),
                                      function(c){c[2]}),
                               function(c){strsplit(c,",")}), 
                        function(vec){as.integer(vec[1])}),
                      var_counts=sapply(
                        sapply(sapply(strsplit(MLS2$RMLS,":"),
                                      function(c){c[2]}),
                               function(c){strsplit(c,",")}), 
                        function(vec){as.integer(vec[2])}))



sample_CN = read.table("RMLS_9_segments.txt", header = T,sep = "\t")

sample_CN  <- sample_CN[,c(1,2,3,11,12)]

names(sample_CN) <- c("chr","start","end","major_cn","minor_cn")

library(GenomicRanges)

library(magrittr)

gsnps = GRanges(seqnames = somatic$chr ,
                ranges = IRanges(somatic$start , somatic$end ),
                strand = "+" )

mcols(gsnps) = somatic


gcnvs = GRanges(seqnames = sample_CN$chr,
                ranges = IRanges(sample_CN$start , sample_CN$end ),
                strand = "+" )


mcols(gcnvs ) = sample_CN

overlaps =findOverlaps(gsnps, gcnvs )

overlaps


merge = cbind( mcols(gsnps[queryHits(overlaps), ]) , mcols(gcnvs[subjectHits(overlaps) ,]))

merge


merge = as.data.frame(merge) %>%
  dplyr::mutate(mutation_id = mutation_id,
                normal_cn = 2,
  ) %>%
  dplyr::select(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn)


merge <- merge[!(merge$major_cn == 0),]

write.table(merge, file = "RMLS_9.tsv", row.names=FALSE, sep="\t",quote = FALSE)

