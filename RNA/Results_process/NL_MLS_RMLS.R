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


Condition <- factor(c(rep("NL",6),rep("MLS",6),rep("RMLS",6)), levels = c("NL","MLS","RMLS"))


colData <- data.frame(row.names = colnames(mycounts),Condition)


colData$Condition <- relevel(colData$Condition, ref = "NL")

colData$Batch <- factor(c(rep("A",3),rep("B",3),rep("A",3), rep("B",3),rep("A",3),rep("B",3)))



expr_count_combat <- ComBat_seq(counts = as.matrix(mycounts), 
                                batch = colData$Batch,
                                group = colData$Condition) 

dds <-DESeqDataSetFromMatrix(expr_count_combat,colData,design = ~ Condition )


dds <- dds[ rowSums(counts(dds)) > 1, ]  


dds <- DESeq(dds)

resultsNames(dds)

res = results(dds,c("Condition", "RMLS","MLS"))


res2 = res[order(res$padj),]

res2


padj.cutoff <- 0.05

cut_lfc <- 1

significant_results <- res2[which(res2$padj < padj.cutoff & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_df <- data.frame(significant_results)

res2 <- data.frame(res2)

sig_MLS_RMLS <- read.table("Sig_DEGs_RMLS_MLS.txt",sep = "\t", header = TRUE)

library(tibble)

significant_results_df <- rownames_to_column(significant_results_df,"Gene")

sig_all <- merge(significant_results_df, sig_MLS_RMLS,"Gene")


library(openxlsx)

write.xlsx(sig_all,"NL_sig_all.xlsx")

##############################################################VENNDIAGRAM

sig_MLS_RMLS_genes <- sig_MLS_RMLS$Gene

sig_NL_MLS_RMLS_genes <- significant_results_df$Gene


dat_gene = list(
  A = sig_MLS_RMLS_genes,
  B = sig_NL_MLS_RMLS_genes)


library(ggVennDiagram)

library(ggpubr)



ggVennDiagram(dat_gene, 
              category.names = c("MLS vs RMLS","NL vs MLS vs RMLS"))+
  scale_fill_gradient(low = "blue", high = "red")+
  ggtitle("Comparison of 2 DEGs sets") + scale_x_continuous(expand = expansion(mult = .6))

library(ggvenn)

names(dat_gene) <- c("MLS vs RMLS","NL vs MLS vs RMLS")


png(filename="final_result/DEGs_venn.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggvenn(dat_gene, columns = c("MLS vs RMLS","NL vs MLS vs RMLS"),
       stroke_size = 0.5,fill_color = c("blue","red"),set_name_size = 5)+
  ggtitle("Comparison of 2 DEGs sets") + theme(plot.title = element_text(size = 15, face = "bold"))

dev.off()

##########################################################
library(DESeq2)

rld <- rlog(dds, blind = FALSE)

rlogMat <- assay(rld)


library(pheatmap)
######################################NO sample names
annotation_col=data.frame(Condition=colData$Condition)

row.names(annotation_col) <- colnames(rlogMat)

annotation_col

################################################

png(filename="final_result/sig_all_NL_DEGs_heatmap.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

pheatmap(rlogMat[sig_all$Gene,],
         annotation_col = annotation_col,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 7,
         scale = "row",
         fontsize_row = 5,
         height = 20,
         annotation_colors = list(Condition=c(NL="green",MLS = "blue",RMLS="red")))

dev.off()
##################################################

intended_gene <- 'CD37'

MLS_gene <- data.frame('values'=rlogMat[intended_gene, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene <- data.frame('values'=rlogMat[intended_gene, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene <- data.frame('values'=rlogMat[intended_gene, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin <- rbind(MLS_gene, RMLS_gene)

for_violin <- rbind(for_violin, NL_gene)

for_violin$Condition <- as.factor(for_violin$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin$Condition  <- factor(for_violin$Condition ,levels = c("NL","MLS","RMLS"))

ggplot(for_violin, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of CD37')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")


###########################################################################

intended_gene2 <- 'SMO'

MLS_gene2 <- data.frame('values'=rlogMat[intended_gene2, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene2 <- data.frame('values'=rlogMat[intended_gene2, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene2 <- data.frame('values'=rlogMat[intended_gene2, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin2 <- rbind(MLS_gene2, RMLS_gene2)

for_violin2 <- rbind(for_violin2, NL_gene2)

for_violin2$Condition <- as.factor(for_violin2$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin2$Condition  <- factor(for_violin2$Condition ,levels = c("NL","MLS","RMLS"))

ggplot(for_violin2, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of SMO')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")


png(filename="final_result/sig_50_DEGs_heatmap.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

intended_gene3 <- 'TP53'

MLS_gene3 <- data.frame('values'=rlogMat[intended_gene3, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene3 <- data.frame('values'=rlogMat[intended_gene3, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene3 <- data.frame('values'=rlogMat[intended_gene3, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin3 <- rbind(MLS_gene3, RMLS_gene3)

for_violin3 <- rbind(for_violin3, NL_gene3)

for_violin3$Condition <- as.factor(for_violin3$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin3$Condition  <- factor(for_violin3$Condition ,levels = c("NL","MLS","RMLS"))



png(filename="final_result/TP53_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(for_violin3, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of TP53')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")

dev.off()



intended_gene4 <- 'IL6'

MLS_gene4 <- data.frame('values'=rlogMat[intended_gene4, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene4 <- data.frame('values'=rlogMat[intended_gene4, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene4 <- data.frame('values'=rlogMat[intended_gene4, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin4 <- rbind(MLS_gene4, RMLS_gene4)

for_violin4 <- rbind(for_violin4, NL_gene4)

for_violin4$Condition <- as.factor(for_violin4$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin4$Condition  <- factor(for_violin4$Condition ,levels = c("NL","MLS","RMLS"))


png(filename="final_result/IL6_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(for_violin4, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of IL6')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")

dev.off()



intended_gene4 <- 'IL6'

MLS_gene4 <- data.frame('values'=rlogMat[intended_gene4, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene4 <- data.frame('values'=rlogMat[intended_gene4, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene4 <- data.frame('values'=rlogMat[intended_gene4, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin4 <- rbind(MLS_gene4, RMLS_gene4)

for_violin4 <- rbind(for_violin4, NL_gene4)

for_violin4$Condition <- as.factor(for_violin4$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin4$Condition  <- factor(for_violin4$Condition ,levels = c("NL","MLS","RMLS"))


png(filename="final_result/IL6_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(for_violin4, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of IL6')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")

dev.off()



intended_gene5 <- 'AKT1'

MLS_gene5 <- data.frame('values'=rlogMat[intended_gene5, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene5 <- data.frame('values'=rlogMat[intended_gene5, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene5 <- data.frame('values'=rlogMat[intended_gene5, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin5 <- rbind(MLS_gene5, RMLS_gene5)

for_violin5 <- rbind(for_violin5, NL_gene5)

for_violin5$Condition <- as.factor(for_violin5$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin5$Condition  <- factor(for_violin5$Condition ,levels = c("NL","MLS","RMLS"))


png(filename="final_result/AKT1_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(for_violin5, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of AKT1')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")

dev.off()


library(dplyr)

intended_gene6 <- 'MICAL3'

MLS_gene6 <- data.frame('values'=rlogMat[intended_gene6, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene6 <- data.frame('values'=rlogMat[intended_gene6, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene6 <- data.frame('values'=rlogMat[intended_gene6, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin6 <- rbind(MLS_gene6, RMLS_gene6)

for_violin6 <- rbind(for_violin6, NL_gene6)

for_violin6$Condition <- as.factor(for_violin6$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin6$Condition  <- factor(for_violin6$Condition ,levels = c("NL","MLS","RMLS"))


png(filename="final_result/MICAL3_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(for_violin6, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of MICAL3')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")

dev.off()

##############################################################################

library(dplyr)

intended_gene7 <- 'IDH2'

MLS_gene7 <- data.frame('values'=rlogMat[intended_gene7, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene7 <- data.frame('values'=rlogMat[intended_gene7, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))

NL_gene7 <- data.frame('values'=rlogMat[intended_gene7, colData$Condition == 'NL'],'Condition'=rep('NL',sum(colData$Condition == 'NL')))


for_violin7 <- rbind(MLS_gene7, RMLS_gene7)

for_violin7 <- rbind(for_violin7, NL_gene7)

for_violin7$Condition <- as.factor(for_violin7$Condition) %>% relevel(ref = 'NL')

library(ggpubr)


for_violin7$Condition  <- factor(for_violin7$Condition ,levels = c("NL","MLS","RMLS"))


png(filename="final_result/IDH2_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(for_violin7, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('NL'='green','MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of IDH2')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "anova")

dev.off()


##############################################################################

library(ggplot2)

library(enrichR)

dbs <- listEnrichrDbs()

dbs <- dbs[order(dbs$libraryName),]

dbs$libraryName

dbs_go <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

dbs_go

dbs_pw <- c("KEGG_2021_Human", "WikiPathways_2023_Human", "Reactome_2022")

dbs_dd <- c("PheWeb_2019", "ClinVar_2019")

dbs_ppi <- c("PPI_Hub_Proteins")


Enriched_go <- enrichr(genes = sig_all$Gene, databases = dbs_go)


Enriched_pw <- enrichr(genes = sig_all$Gene, databases = dbs_pw)


Enriched_ppi <- enrichr(genes = sig_all$Gene, databases = dbs_ppi)


plotEnrich(Enriched_ppi[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value")

Enriched_ppi$PPI_Hub_Proteins

#####################################################################################
library(tinyarray)
library(tidyverse)

library(rtracklayer)


Gle <- read.table("MLS1.txt",skip = 1,sep="\t",header = T)

Gle  <- Gle[,c(1,6)]

le = Gle[match(rownames(expr_count_combat),Gle$Geneid),"Length"]

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

tpms <- apply(expr_count_combat,2,countToTpm,le)


exp = as.data.frame(tpms)

exp = rownames_to_column(exp)

write.table(exp,file = "exp_NL_MLS_RMLS_tpm.txt",row.names = F,quote = F,sep = "\t")

source("CIBERSORT.R")

TME = CIBERSORT("LM22.txt", 
                "exp_NL_MLS_RMLS_tpm.txt", 
                perm = 100, 
                QN = T)

cibersort_res <- TME[,-(23:25)]

library(pheatmap)

k <- apply(cibersort_res,2,function(x) {sum(x == 0) < nrow(TME)/2})

table(k)

cibersort_res2 <- as.data.frame(t(cibersort_res[,k]))


ann_colr <- list(C = c(MLS = "green",RMLS = "red"))


png(filename="final_result/CIBERSORT_NL_heatmap_re.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

pheatmap(cibersort_res2,scale = "row",
         show_colnames = T,
         cluster_cols = F,
         annotation_col = colData,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         annotation_colors = list(Batch=c(A = "#D95F02",B="#7570B3"),
           Condition=c(NL ="#1B9E77",MLS ="#2166AC",RMLS="#C11C0F")))

dev.off()

dat <- cibersort_res%>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

dat$Condition = ifelse(as.character(str_sub(dat$Sample)) %in% c("NL3","NL4","NL5","NL6","NL7","NL8"),"NL",ifelse(as.character(str_sub(dat$Sample)) %in% c("MLS3","MLS4","MLS5","MLS6","MLS7","MLS8"),"MLS","RMLS"))



library(ggpubr)

library(dplyr)

library(tidyr)

library(tibble)

library(ggplot2)

library(stringr)

library(ggsignif)


ggplot(dat,aes(Cell_type,Proportion,fill = Condition)) + 
  geom_boxplot(color = "black") + 
  theme_bw() + 
  labs(x = "Immune cell type", y = "Proportion of Immune signature") +
  theme(legend.position = "right") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = c('NL'= 'green','MLS'='blue', 'RMLS'='red'))+ stat_compare_means(aes(group = Condition,label = ..p.signif..),method = "anova")

