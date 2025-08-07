library(DESeq2)

library(tidyverse)

library(sva)

library(tibble)

NL3 <- read.table("NL3.txt", sep = "\t", header = T)
NL4 <- read.table("NL4.txt", sep = "\t", header = T)
NL5 <- read.table("NL5.txt", sep = "\t", header = T)
NL6 <- read.table("15190-NL.txt",sep = "\t",header = T)
NL7 <- read.table("57257-NL.txt",sep = "\t",header = T)
NL8 <- read.table("70615-NL.txt",sep = "\t",header = T)



colnames(NL3)[7]<-("NL3")
colnames(NL4)[7]<-("NL4")
colnames(NL5)[7]<-("NL5")
colnames(NL6)[7]<-("NL6")
colnames(NL7)[7]<-("NL7")
colnames(NL8)[7]<-("NL8")


NL3 <- NL3[,-(2:6)]
NL4 <- NL4[,-(2:6)]
NL5 <- NL5[,-(2:6)]
NL6 <- NL6[,-(2:6)]
NL7 <- NL7[,-(2:6)]
NL8 <- NL8[,-(2:6)]

n_raw_count <- merge(NL3,NL4,by="Geneid")
n_raw_count2 <- merge(NL5,NL6,by="Geneid")
n_raw_count3 <- merge(NL7,NL8,by="Geneid")


n_raw_count <- merge(n_raw_count,n_raw_count3,by="Geneid")

mycounts <- read.csv("7RMLS_7MLS_readout2.csv",sep = ",")

mycounts =  merge(n_raw_count,mycounts, by = "Geneid")


mycounts<-mycounts[,-c(14,21)]

rownames(mycounts)<-mycounts[,1]

mycounts<-mycounts[,-1]

write.table(mycounts,"NL_MLS_RMLS_raw_counts.txt",sep = "\t",quote = FALSE)


mycounts <- read.table("NL_MLS_RMLS_raw_counts.txt",sep = "\t",header = TRUE)

mycounts <- mycounts[,-c(7:12)]



Condition <- factor(c(rep("NL",6),rep("RMLS",6)), levels = c("NL","RMLS"))


colData <- data.frame(row.names = colnames(mycounts),Condition)


colData$Condition <- relevel(colData$Condition, ref = "NL")

colData$Batch <- factor(c(rep("A",3),rep("B",3),rep("A",3),rep("B",3)))

expr_count_combat <- ComBat_seq(counts = as.matrix(mycounts), 
                                batch = colData$Batch ) 

dds <-DESeqDataSetFromMatrix(expr_count_combat,colData,design = ~ Condition )


dds <- dds[ rowSums(counts(dds)) > 1, ]  


dds <- DESeq(dds)

resultsNames(dds)

res = results(dds,c("Condition", "RMLS","NL"))

res2 = res[order(res$padj),]

res2


padj.cutoff <- 0.05

cut_lfc <- 1

significant_results <- res2[which(res2$padj < padj.cutoff & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_df <- data.frame(significant_results)

res2 <- data.frame(res2)

sig_MLS_RMLS <- read.table("Sig_DEGs_RMLS_MLS.txt",sep = "\t", header = TRUE)

sig_MLS_RMLS 


library(tibble)

significant_results_df <- rownames_to_column(significant_results_df,"Gene")


write.table(significant_results_df,"Sig_DEGs_RMLS_NL.txt",sep = "\t",quote = FALSE,row.names = FALSE)

sig_all <- merge(significant_results_df, sig_MLS_RMLS,"Gene")

colnames(sig_all) <- c("Gene","baseMean_vsNL","log2FoldChange_vsNL",
                       "lfcSE_vsNL","stat_vsNL","pvalue_vsNL",
                       "padj_vsNL",
                       "baseMean_vsMLS","log2FoldChange_vsMLS",
                       "lfcSE_vsMLS","stat_vsMLS","pvalue_vsMLS",
                       "padj_vsMLS")

library(openxlsx)

write.xlsx(sig_all,"NL_RMLS_MLS_sig_all.xlsx")



sig_MLS_RMLS_genes <- sig_MLS_RMLS$Gene

sig_NL_MLS_RMLS_genes <- significant_results_df$Gene


dat_gene = list(
  A = sig_MLS_RMLS_genes,
  B = sig_NL_MLS_RMLS_genes)


library(ggVennDiagram)

library(ggpubr)



ggVennDiagram(dat_gene, 
              category.names = c("MLS vs RMLS","NL vs RMLS"))+
  scale_fill_gradient(low = "blue", high = "red")+
  ggtitle("Comparison of 2 DEGs sets") + scale_x_continuous(expand = expansion(mult = .6))

library(ggvenn)

names(dat_gene) <- c("MLS vs RMLS","NL vs RMLS")


png(filename="final_result/DEGs_venns_re.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggvenn(dat_gene, columns = c("MLS vs RMLS","NL vs RMLS"),
       stroke_size = 0.5,fill_color = c("blue","red"),set_name_size = 5)+
  ggtitle("Comparison of 2 DEGs sets") + theme(plot.title = element_text(size = 15, face = "bold"))

dev.off()
