library(DESeq2)

library(tidyverse)

library(sva)

library(tibble)

mycounts<-read.csv("7RMLS_7MLS_readout2.csv")

mycounts<-mycounts[,-c(8,15)]


rownames(mycounts)<-mycounts[,1]

mycounts<-mycounts[,-1]

Condition <- factor(c(rep("MLS",6),rep("RMLS",6)), levels = c("MLS","RMLS"))


colData <- data.frame(row.names = colnames(mycounts),Condition)

colData$Condition <- relevel(colData$Condition, ref = "MLS")

colData$Batch <- factor(c(rep("A",3),rep("B",3),rep("A",3), rep("B",3)))


expr_count_combat <- ComBat_seq(counts = as.matrix(mycounts), 
                                batch = colData$Batch  ) 


dds <-DESeqDataSetFromMatrix(expr_count_combat,colData,design = ~ Condition)


dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq(dds)

res = results(dds,c("Condition", "RMLS","MLS"))

res2 = res[order(res$padj),]

res2

rld <- rlog(dds, blind = FALSE)

rlogMat <- assay(rld)


padj.cutoff <- 0.05

cut_lfc <- 1

significant_results <- res2[which(res2$padj < padj.cutoff & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results = significant_results[order(significant_results$padj),]

significant_results_df <- rownames_to_column(data.frame(significant_results),"Gene")

write.table(significant_results_df,"Sig_DEGs_RMLS_MLS.txt",sep = "\t",quote = FALSE,row.names = FALSE)

significant_results_High <- res2[which(res2$padj < padj.cutoff & (res2$log2FoldChange>cut_lfc)),]

significant_results_High_genes <- rownames(significant_results_High)

significant_results_Low <- res2[which(res2$padj < padj.cutoff & (res2$log2FoldChange<(-cut_lfc) )),]


significant_results_Low_genes <- rownames(significant_results_Low)



significant_genes_50 <- rownames(significant_results[1:50,])


significant_genes_100 <- rownames(significant_results[1:100,])


significant_results_High_25 <- significant_results_High[1:25,]

significant_results_Low_25 <- significant_results_Low[1:25,]

significant_results_HL_50 <- rbind(significant_results_High_25,significant_results_Low_25)

significant_results_HL_50_genes <- rownames(significant_results_HL_50)


library(pheatmap)
######################################annotation sample names
annotation_col=data.frame(Condition=colData$Condition)

row.names(annotation_col) <- colnames(rlogMat)

annotation_col

################################################

library(pheatmap)

png(filename="final_result/sig_50_DEGs_heatmap.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

pheatmap(rlogMat[significant_genes_50,],
         annotation_col = annotation_col,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 7,
         scale = "row",
         fontsize_row = 5,
         height = 20,
         annotation_colors = list(Condition=c(MLS = "blue",RMLS="red")))

dev.off()


pheatmap(rlogMat[c("HES5","THBS2","ZNF460"),],
         annotation_col = annotation_col,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 7,
         scale = "row",
         fontsize_row = 5,
         height = 20,
         annotation_colors = list(Condition=c(MLS = "blue",RMLS="red")))


topT <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(topT)[1] <- "Gene"

cut_pvalue <- 0.05


png(filename="final_result/sig_50_DEGs_volcano.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

par(mar=c(5,5,5,2))
# Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="RMLS vs MLS Volcano plot", col='grey',  xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Adj~p~value), lwd = 2.0,cex = 3,cex.lab = 1.5,cex.axis = 1.5,cex.main=2,ylim = c(0,10)))

with(topT,legend("right",c("MLS","Not sig","RMLS"),col = c('blue',"grey",'red'),pch = 20, title = "Condition"))


with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc),points(log2FoldChange, -log10(padj), pch=20, col='red', cex=2))

with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=2))

## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=2.0)

abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)

abline(v=cut_lfc, col='black', lty=4, lwd=2.0)

abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)


library(calibrate)

topTT <- column_to_rownames(topT,"Gene")


with(topTT[significant_results_HL_50_genes, ], {
  text(log2FoldChange, -log10(padj), significant_results_HL_50_genes, pos=2,font = 2,cex=.5)
})

dev.off()
###############################################################################################################################
with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), textxy(log2FoldChange, -log10(padj), labs=Gene,font = 2, cex=1))

with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), textxy(log2FoldChange, -log10(padj), labs=Gene,font = 2, cex=1))
##################################################################################################################################
###################################################################################FOR_GSEA_Prelanked
for_gsea_pre = data.frame(res2)

for_gsea_pre = rownames_to_column(for_gsea_pre,"Gene")

for_gsea_pre$fcsign = sign(for_gsea_pre$log2FoldChange)

for_gsea_pre$logP = -log10(for_gsea_pre$pvalue)

for_gsea_pre$metric = for_gsea_pre$logP/for_gsea_pre$fcsign

for_gsea = for_gsea_pre[,c("Gene","metric")]

for_gsea = na.omit(for_gsea)

write.table(for_gsea,file="DE_sig_genes.rnk",quote=F,sep="\t",row.names=F,col.names = F)


###################################################################FOR_CIBERSORT
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

write.table(exp,file = "exp_tpm.txt",row.names = F,quote = F,sep = "\t")

source("CIBERSORT.R")

TME = CIBERSORT("LM22.txt", 
                "exp_DESeq.txt", 
                perm = 100, 
                QN = T)

cibersort_res <- TME[,-(23:25)]

library(pheatmap)

k <- apply(cibersort_res,2,function(x) {sum(x == 0) < nrow(TME)/2})

table(k)

cibersort_res2 <- as.data.frame(t(cibersort_res[,k]))


ann_colr <- list(C = c(MLS = "green",RMLS = "red"))

pheatmap(cibersort_res2,scale = "row",
         show_colnames = T,
         cluster_cols = F,
         annotation_col = colData,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         annotation_colors = list(Condition=c(MLS = "blue",RMLS="red")))

###################################################################################################FOR_ESTIMATE


normalized_counts <- counts(dds, normalized=TRUE)

Res_count_data = data.frame(normalized_counts)

Res_count_data = rownames_to_column(Res_count_data,"Gene")


write.table(Res_count_data ,file = "exp_DESeq.txt",row.names = F,quote = F,sep = "\t")


library(estimate)

library(ggplot2)

library(ggsignif)


filterCommonGenes(input.f = "exp_DESeq.txt",
                  output.f = "MLS_RMLS_genes.gct",
                  id= "GeneSymbol")

d = estimate::SI_geneset


estimateScore(input.ds="MLS_RMLS_genes.gct", output.ds="MLS_RMLS_estimate_score.gct", platform="affymetrix")

estimate_score <- read.table("MLS_RMLS_estimate_score.gct", skip = 2, header = TRUE)


library(dplyr)

estimate_score_df =  rownames_to_column(estimate_score,"Sample")

write.csv(estimate_score_df,"ESTIMATE_score_deseq2.csv",row.names = FALSE)


rownames(estimate_score) = estimate_score[,1]

estimate_score = t(estimate_score[,3:ncol(estimate_score)])

estimate_score = data.frame(estimate_score)

estimate_score$Condition = c(rep("MLS",6),rep("RMLS",6))

estimate_score$Condition<- as.factor(estimate_score$Condition) %>% relevel(ref = 'MLS')


png(filename="final_result/stromalscore.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

ggplot(estimate_score, aes(x=Condition, y=StromalScore, fill=Condition))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c('MLS'='blue', 'RMLS'='red'))+
  geom_signif(comparisons = list(c("MLS","RMLS")),
              map_signif_level = FALSE,test = "t.test")+
  labs(x= 'Condition', y= 'StromalScore')

dev.off()


png(filename="final_result/immunescore.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggplot(estimate_score, aes(x=Condition, y=ImmuneScore, fill=Condition))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c('MLS'='blue', 'RMLS'='red'))+
  geom_signif(comparisons = list(c("MLS","RMLS")),
              map_signif_level = FALSE,test = "t.test")+
  labs(x= 'Condition', y= 'ImmuneScore')

dev.off()


png(filename="final_result/estimatescore.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

ggplot(estimate_score, aes(x=Condition, y=ESTIMATEScore, fill=Condition))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c('MLS'='blue', 'RMLS'='red'))+
  geom_signif(comparisons = list(c("MLS","RMLS")),
              map_signif_level = FALSE,test = "t.test")+
  labs(x= 'Condition', y= 'ESTIMATEScore')

dev.off()


png(filename="final_result/tumorpurity.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

ggplot(estimate_score, aes(x=Condition, y=TumorPurity, fill=Condition))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = TRUE)+
  scale_fill_manual(values = c('MLS'='blue', 'RMLS'='red'))+
  geom_signif(comparisons = list(c("MLS","RMLS")),
              map_signif_level = FALSE,test = "t.test")+
  labs(x= 'Condition', y= 'TumorPurity')

dev.off()

################################################################################################

Sonic_Hedgehog_signaling  = c("HHAT","HHIP","DHH","IHH",
                              "HES1","HES2","HEY1","CCND1")

#######################################################################################
library(ggplot2)

library(enrichR)

dbs <- listEnrichrDbs()

dbs <- dbs[order(dbs$libraryName),]

dbs$libraryName

dbs_go <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

dbs_go

dbs_pw <- c("KEGG_2021_Human", "WikiPathways_2021_Human", "Reactome_2022")

dbs_dd <- c("PheWeb_2019", "ClinVar_2019")


dbs_hn <- c("PheWeb_2019", "ClinVar_2019")


dbs_ppi <- c("PPI_Hub_Proteins")


Enriched_high_go <- enrichr(genes = significant_results_High_genes, databases = dbs_go)



Enriched_low_go <- enrichr(genes = significant_results_Low_genes, databases = dbs_go)



plotEnrich(Enriched_high_go[[3]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value")


Enriched_high_pw <- enrichr(genes = significant_results_High_genes, databases = dbs_pw)



Enriched_low_pw <- enrichr(genes = significant_results_Low_genes, databases = dbs_pw)


plotEnrich(Enriched_high_pw[[2]], showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value")


plotEnrich(Enriched_low_pw[[1]], showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value")


Enriched_high_ppi <- enrichr(genes = significant_results_High_genes, databases = dbs_ppi)



plotEnrich(Enriched_high_ppi[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value")



Enriched_low_ppi <- enrichr(genes = significant_results_Low_genes, databases = dbs_ppi)



plotEnrich(Enriched_low_ppi[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value")

library(Rtoolbox)



replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_CELL_CYCLE",class.name = "Round cell MLS")


replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_HEDGEHOG_SIGNALING_PATHWAY",class.name = "Round cell MLS")


replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_NOTCH_SIGNALING_PATHWAY",class.name = "Round cell MLS")


replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_CITRATE_CYCLE_TCA_CYCLE",class.name = "Round cell MLS")

replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_GLUTATHIONE_METABOLISM",class.name = "Round cell MLS")

replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_ECM_RECEPTOR_INTERACTION",class.name = "Round cell MLS")


replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_GAP_JUNCTION",class.name = "Round cell MLS")


replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725528080990(KEGG)/","KEGG_CELL_ADHESION_MOLECULES_CAMS",class.name = "Round cell MLS")



replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725526490749(GO_BP)/","GOBP_KERATINOCYTE_DIFFERENTIATION",class.name = "Round cell MLS")



replotGSEA("RMLS_MLS_DESEQ.GseaPreranked.1725526490749(GO_BP)/","GOBP_EPIDERMAL_CELL_DIFFERENTIATION",class.name = "Round cell MLS")

##################################################################################################

library(dplyr)

intended_gene <- 'AKT1'

MLS_gene <- data.frame('values'=rlogMat[intended_gene, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene <- data.frame('values'=rlogMat[intended_gene, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))


for_violin <- rbind(MLS_gene, RMLS_gene)


for_violin$Condition <- as.factor(for_violin$Condition) %>% relevel(ref = 'MLS')

library(ggpubr)


for_violin$Condition  <- factor(for_violin$Condition ,levels = c("MLS","RMLS"))


png(filename="final_result/AKT1_2DGs_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)



ggplot(for_violin, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of AKT1')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "t.test")

dev.off()

##################################################################################################

library(dplyr)

intended_gene2 <- 'PATZ1'

MLS_gene2 <- data.frame('values'=rlogMat[intended_gene2, colData$Condition == 'MLS'],'Condition'=rep('MLS',sum(colData$Condition == 'MLS')))

RMLS_gene2 <- data.frame('values'=rlogMat[intended_gene2, colData$Condition == 'RMLS'],'Condition'=rep('RMLS',sum(colData$Condition == 'RMLS')))


for_violin2 <- rbind(MLS_gene2, RMLS_gene2)


for_violin2$Condition <- as.factor(for_violin2$Condition) %>% relevel(ref = 'MLS')

library(ggpubr)


for_violin2$Condition  <- factor(for_violin2$Condition ,levels = c("MLS","RMLS"))


png(filename="final_result/PATZ1_2DGs_expression_violin.png",width=1000,height=800,unit="px",bg="transparent",res = 150)



ggplot(for_violin2, aes(x=Condition, y=values, fill=Condition))+
  geom_violin(position = position_dodge(width = 1),scale = 'width')+
  scale_fill_manual(values = c('MLS'='blue', 'RMLS'='red'))+
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+
  labs(x= 'Condition', y= 'normalized expression counts of PATZ1')+theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size = 1,fill = 'white'),
        panel.grid = element_blank()) + stat_compare_means(method = "t.test")

dev.off()

###############################################################################CIBERSORT

png(filename="final_result/CIBERSORT_NL_heatmap_re.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

pheatmap(cibersort_res2,scale = "row",
         show_colnames = T,
         cluster_cols = F,
         annotation_col = colData,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         annotation_colors = list(Batch=c(A = "#D95F02",B="#7570B3"),
                                  Condition=c(NL ="#1B9E77",MLS ="#2166AC",RMLS="#C11C0F")))

dev.off()
