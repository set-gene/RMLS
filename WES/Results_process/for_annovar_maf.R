
MLS1 <- read.table("MLS_1.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS1$Tumor_Sample_Barcode = "MLS_1"

MLS2 <-  read.table("MLS_2.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS2$Tumor_Sample_Barcode = "MLS_2"

MLS3<-  read.table("MLS_3.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS3$Tumor_Sample_Barcode = "MLS_3"

MLS4 <- read.table("MLS_4.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS4$Tumor_Sample_Barcode = "MLS_4"

MLS5 <- read.table("MLS_5.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS5$Tumor_Sample_Barcode = "MLS_5"

MLS6 <-  read.table("MLS_6.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS6$Tumor_Sample_Barcode = "MLS_6"

MLS7 <-  read.table("MLS_7.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS7$Tumor_Sample_Barcode = "MLS_7"


MLS8 <-  read.table("MLS_8.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS8$Tumor_Sample_Barcode = "MLS_8"

MLS9 <-  read.table("MLS_9.hg38_multianno.txt",sep = "\t",header = TRUE)

MLS9$Tumor_Sample_Barcode = "MLS_9"


library(dplyr)

MLS <-  bind_rows(MLS1,MLS2,MLS3,MLS4,MLS5,MLS6,MLS7,MLS8,MLS9) 


write.table(MLS,"merge_MLS_annovar.txt",sep = "\t",quote = FALSE, row.names = FALSE)

library(maftools)
###########################################################vaf
MLS.var.annovar2.maf = annovarToMaf(annovar = "merge_MLS_annovar.txt", 
                                   Center = 'NA', 
                                   refBuild = 'hg38', 
                                   tsbCol = 'Tumor_Sample_Barcode', 
                                   table = 'refGene',
                                   sep = "\t")
library(ICAMS)

MLS.var.annovar2.maf$Variant_Classification

df2 <-MLS.var.annovar2.maf[,c(69:79)]

colnames(df2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","MLS","NL")

df_3_vaf <- GetMutectVAF(df2,tumor.col.name = "MLS")


MLS.var.annovar2.maf <-MLS.var.annovar2.maf[,-c(69:79)]


MLS.var.annovar2.maf <- cbind(MLS.var.annovar2.maf, df_3_vaf)


write.table(MLS.var.annovar2.maf, file="MLS.WES.merge.maf",quote=F,sep = "\t",row.names = F)

############################################################################################################
RMLS1 <- read.table("RMLS_1.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS1$Tumor_Sample_Barcode = "RMLS_1"

RMLS2 <-  read.table("RMLS_2.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS2$Tumor_Sample_Barcode = "RMLS_2"

RMLS3<-  read.table("RMLS_3.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS3$Tumor_Sample_Barcode = "RMLS_3"

RMLS4 <- read.table("RMLS_4.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS4$Tumor_Sample_Barcode = "RMLS_4"

RMLS5 <- read.table("RMLS_5.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS5$Tumor_Sample_Barcode = "RMLS_5"

RMLS6 <-  read.table("RMLS_6.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS6$Tumor_Sample_Barcode = "RMLS_6"

RMLS7 <-  read.table("RMLS_7.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS7$Tumor_Sample_Barcode = "RMLS_7"

RMLS8 <-  read.table("RMLS_8.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS8$Tumor_Sample_Barcode = "RMLS_8"

RMLS9 <-  read.table("RMLS_9.hg38_multianno.txt",sep = "\t",header = TRUE)

RMLS9$Tumor_Sample_Barcode = "RMLS_9"

library(dplyr)

RMLS <-  bind_rows(RMLS1,RMLS2,RMLS3,RMLS4,RMLS5,RMLS6,RMLS7,RMLS8,RMLS9) 


write.table(RMLS,"merge_RMLS_annovar.txt",sep = "\t",quote = FALSE, row.names = FALSE)

library(maftools)
###########################################################vaf
RMLS.var.annovar2.maf = annovarToMaf(annovar = "merge_RMLS_annovar.txt", 
                                    Center = 'NA', 
                                    refBuild = 'hg38', 
                                    tsbCol = 'Tumor_Sample_Barcode', 
                                    table = 'refGene',
                                    sep = "\t")
library(ICAMS)

RMLS.var.annovar2.maf$Variant_Classification

Rdf2 <-RMLS.var.annovar2.maf[,c(69:79)]

colnames(Rdf2) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NL","RMLS")

Rdf_3_vaf <- GetMutectVAF(Rdf2,tumor.col.name = "RMLS")


RMLS.var.annovar2.maf <-RMLS.var.annovar2.maf[,-c(69:79)]


RMLS.var.annovar2.maf <- cbind(RMLS.var.annovar2.maf, Rdf_3_vaf)


write.table(RMLS.var.annovar2.maf, file="RMLS.WES.merge.maf",quote=F,sep = "\t",row.names = F)
