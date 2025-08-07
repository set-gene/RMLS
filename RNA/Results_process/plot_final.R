#volcano plot 

library(ggplot2)

library(ggpubr)

library(ggrepel)

library(openxlsx)
library(RColorBrewer)

brewer.pal( name = "RdBu",n = 11)


sig_all <- read.xlsx("NL_RMLS_MLS_sig_all.xlsx")


sig_all_high <- subset(sig_all, sig_all$log2FoldChange_vsNL > cut_lfc & sig_all$log2FoldChange_vsMLS > cut_lfc)

sig_all_high_genes <- sig_all_high$Gene

sig_all_low <- subset(sig_all, sig_all$log2FoldChange_vsNL <(-cut_lfc) & sig_all$log2FoldChange_vsMLS <(-cut_lfc))

sig_all_low_genes <- sig_all_low$Gene

topT$DEGs <- ifelse(topT$Gene %in% sig_all_high_genes,"UP",ifelse(topT$Gene %in% sig_all_low_genes,"DOWN","NOT"))

topT$sig_gene <- ifelse(topT$Gene %in% c("CD37"), topT$Gene,NA)


png(filename="final_result/VOLCANO_final.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(data = topT, aes(x=log2FoldChange,y=-log10(padj),col = DEGs)) +
  geom_point(alpha=0.4 ,shape = 16)+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept = -log10(padj.cutoff) ,linetype=4) +
  geom_vline(xintercept = c(-1,0,1) ,linetype=4) +
  scale_color_manual(name = "", values = c("#C11C0F","#2166AC", "gray50"), limits = c("UP", "DOWN", "NOT")) +
  geom_text_repel(aes(label=sig_gene), fontface="bold", color="black",point.padding = 0.2,nudge_x = .15,
                  nudge_y = .5,
                  segment.linetype = 6,
                  segment.curvature = -1e-20,
                  arrow = arrow(length = unit(0.015, "npc"))) + theme(axis.text.x = element_text(family = "Arial",size = 9,color = "black",vjust = 0.5, hjust=1,face = "bold"),                  
        axis.text.y = element_text(family  = "Arial",size = 10,color = "black",face = "bold"),                   # face the y axis title/label
        axis.title.y = element_text(family = "Arial",hjust = 0.4,face = "bold",size = 12),                   # face the y axis title/label
        axis.title.x = element_text(family = "Arial",hjust = 0.5,face = "bold",size = 12)) + ylab("-Log10(Adj.P.value)") + xlab("Log2FoldChange")

dev.off()

#VENNDIAGRAM 

sig_MLS_RMLS <- read.table("Sig_DEGs_RMLS_MLS.txt",sep = "\t", header = TRUE,stringsAsFactors = FALSE)


sig_NL_RMLS <- read.table("Sig_DEGs_RMLS_NL.txt",sep = "\t", header = TRUE,stringsAsFactors = FALSE)


sig_MLS_RMLS_genes <- sig_MLS_RMLS$Gene

sig_NL_RMLS_genes <- sig_NL_RMLS$Gene



dat_gene = list(A = sig_MLS_RMLS_genes, B = sig_NL_RMLS_genes)


library(ggVennDiagram)

library(ggpubr)

library(ggplot2)

library(officer)
library(rvg)

pvenn <- ggVennDiagram(dat_gene, 
              label_alpha = 0,
              category.names = c("MLS vs RMLS","NL vs RMLS"),
              set_color = c("#D55E00","#CC79A7"),) +
  scale_x_continuous(expand = expansion(mult = .6)) + scale_fill_distiller(palette = "Reds", direction = 1)
#####################################
read_pptx() %>% 
  print(target = "./example.pptx")

read_pptx( "./example.pptx")


editable_graph <- rvg::dml(ggobj = pvenn)

read_pptx( "./example.pptx") %>%
  add_slide("Title Slide","Office Theme") %>%
  ph_with(editable_graph,
          location = ph_location(left=0, top=0,width = 5, height = 5, bg="transparent")) %>%
  print(target = "./example.pptx")

###################################
#########box&violin

all_tpm <- read.table("exp_NL_MLS_RMLS_tpm.txt",sep = "\t",stringsAsFactors = FALSE, header = TRUE)

library(tibble)

all_tpm <- column_to_rownames(all_tpm,"rowname")

all_log_tpm <- log2(all_tpm +1)

all_tpm_t <- as.data.frame(t(all_log_tpm))

select_tpm <- all_tpm_t[,c("AKT1","PATZ1","SMO","TP53","IL6","MICAL3","TNF")]


Condition2 <- factor(c(rep("NL",6),rep("MLS",6),rep("RMLS",6)), levels = c("NL","MLS","RMLS"))


colData2 <- data.frame(row.names = colnames(all_tpm),Condition2)


colData2$Condition2 <- relevel(colData2$Condition2, ref = "NL")

colData2$Batch <- factor(c(rep("A",3),rep("B",3),rep("A",3), rep("B",3),rep("A",3),rep("B",3)))

select_tpm$Condition <- colData2$Condition2

#data_summary <- function(x) {
 # m <- mean(x)
 # ymin <- m-sd(x)
 # ymax <- m+sd(x)
  #return(c(y=m,ymin=ymin,ymax=ymax))
#}


png(filename="final_result/AKT1_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=AKT1, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of AKT1") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()



png(filename="final_result/SMO_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=SMO, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of SMO") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()


png(filename="final_result/TP53_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=TP53, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of TP53") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()


png(filename="final_result/IL6_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=IL6, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of IL6") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()



png(filename="final_result/TNF_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=TNF, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of TNF") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()



png(filename="final_result/PATZ1_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=PATZ1, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of PATZ1") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()




png(filename="final_result/MICAL3_EXP_violin.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(select_tpm, aes(x=Condition, y=MICAL3, fill=Condition)) + 
  geom_violin()+
  scale_y_continuous(name = "Normalized Expression Counts of MICAL3") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#1B9E77","#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "anova") 

dev.off()
#AKT1, PATZ1, SMO,TP53  ,IL6
#PATHway cell cycle Notch hedgehog DNA repair fattyacid ERBB 
#ECM 

#GO_BP_GSEA_up <- read.table("gsea_report_for_na_pos_1725526490749_GOBP.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE)

#GO_BP_GSEA_up <- GO_BP_GSEA_up[,c(1,6,8)]

#GO_BP_GSEA_up_sig <- GO_BP_GSEA_up[which(GO_BP_GSEA_up$FDR.q.val < 0.25 & (GO_BP_GSEA_up$NES> 1)),]

#GO_BP_GSEA_up_sig_select <- subset(GO_BP_GSEA_up_sig, GO_BP_GSEA_up_sig %in% "CELL_CYCLE"| GO_BP_GSEA_up_sig %in% "NOTCH"|GO_BP_GSEA_up_sig %in%  "HEDGEHOG"| GO_BP_GSEA_up_sig %in% "DNA_REPAIR"| GO_BP_GSEA_up_sig %in% "DOUBLE_STRAND_BREAK_REPAIR")


KEGG_GSEA_up <- read.table("gsea_report_for_na_pos_1725528080990_KEGG.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE)


KEGG_GSEA_up <- KEGG_GSEA_up[,c(1,6,8)]

KEGG_GSEA_up_sig <- KEGG_GSEA_up[which(KEGG_GSEA_up$FDR.q.val < 0.25 & (KEGG_GSEA_up$NES> 1)),]

library(dplyr)

KEGG_GSEA_up_sig_select <- KEGG_GSEA_up_sig %>% filter(grepl("CELL_CYCLE",NAME)| grepl("NOTCH",NAME) | grepl("HEDGEHOG",NAME) | grepl("REPAIR",NAME) | grepl("MTOR",NAME)| grepl("P53",NAME) |  grepl("REPAIR",NAME) | grepl("ACUTE_MYELOID_LEUKEMIA",NAME) | grepl("FATTY_ACID_METABOLISM" ,NAME)) 

theTable <- within(KEGG_GSEA_up_sig_select, 
                   NAME <- factor(NAME,levels=KEGG_GSEA_up_sig_select[order(KEGG_GSEA_up_sig_select$NES,decreasing=T),]$NAME))


mid <- mean(theTable$FDR.q.val)

library(ggthemes)


png(filename="final_result/KEGG_UP_REG.png",width=800,height=600,unit="px",bg="transparent",res = 100)



ggplot(theTable, aes(y=NES, x=NAME)) +
  geom_bar(stat="identity",aes(fill = FDR.q.val))+
  scale_fill_gradient(low ="#CC79A7", high = "#D55E00")+
  ylab("Normalized Enrichment Score")+
  xlab("")+
  geom_hline(yintercept = 1.0,linetype='dashed')+
  theme_few()+
  scale_y_continuous(limits = c(0, 2.5))+
  theme(axis.text.x = element_text(face = "bold",angle = 75,vjust = 1,hjust = 1,size = 8,family = "arial",colour = "black"),                      # adjusting the position
        axis.title.y = element_text(face = "bold",size = 8,family = "arial"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 8,family = "arial",colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 10))

dev.off()

KEGG_GSEA_dn <- read.table("gsea_report_for_na_neg_1725528080990_KEGG.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE)


KEGG_GSEA_dn <- KEGG_GSEA_dn[,c(1,6,8)]

KEGG_GSEA_dn_sig <- KEGG_GSEA_dn[which(KEGG_GSEA_dn$FDR.q.val < 0.25 & KEGG_GSEA_dn$NES < (-1)),]


library(dplyr)

KEGG_GSEA_dn_sig_select <- KEGG_GSEA_dn_sig %>% filter(grepl("CYTOKINE",NAME)| grepl("JAK_STAT",NAME) | grepl("P450",NAME) | grepl("T_CELL",NAME) | grepl("CHEMOKINE",NAME)| grepl("ECM",NAME) |  grepl("GAP",NAME) | grepl("ADHESION",NAME) | grepl("TGF" ,NAME)) 


theTable2 <- within(KEGG_GSEA_dn_sig_select, 
                   NAME <- factor(NAME,levels=KEGG_GSEA_dn_sig_select[order(KEGG_GSEA_dn_sig_select$NES,decreasing=T),]$NAME))


png(filename="final_result/KEGG_DN_REG.png",width=800,height=600,unit="px",bg="transparent",res = 100)


ggplot(theTable2, aes(y=NES, x=NAME)) +
  geom_bar(stat="identity",aes(fill = FDR.q.val))+
  scale_fill_gradient(low ="#CC79A7", high = "#D55E00")+
  ylab("Normalized Enrichment Score")+
  xlab("")+
  geom_hline(yintercept = -1.0,linetype='dashed')+
  scale_y_continuous(limits = c(-2.5, 0))+
  theme_few()+
  theme(axis.text.x = element_text(face = "bold",angle = 75,vjust = 1,hjust = 1,size = 8,family = "arial",colour = "black"),                      # adjusting the position
        axis.title.y = element_text(face = "bold",size = 8,family = "arial"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 8,family = "arial",colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 10))

dev.off()

theTable3 <- rbind(theTable,theTable2)

library(ggthemes)
library(ggplot2)
library(ggpubr)

png(filename="final_result/KEGG_all_REG.png",width=800,height=600,unit="px",bg="transparent",res = 100)

ggplot(theTable3, aes(y=NES, x=NAME)) +
  geom_bar(stat="identity",aes(fill = FDR.q.val))+
  scale_fill_gradient(low ="#C11C0F", high = "#2166AC")+
  ylab("Normalized Enrichment Score")+
  xlab("")+
  geom_hline(yintercept = 1.0,linetype='dashed')+
  geom_hline(yintercept = -1.0,linetype='dashed')+
  scale_y_continuous(limits = c(-2.5, 2.5))+
  theme_few()+
  theme(axis.text.x = element_text(face = "bold",angle = 75,vjust = 1,hjust = 1,size = 8,family = "arial",colour = "black"),                      # adjusting the position
        axis.title.y = element_text(face = "bold",size = 8,family = "arial"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 8,family = "arial",colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 10))

dev.off()

GO_GSEA_dn <- read.table("gsea_report_for_na_neg_1725526490749_GO_BP.tsv",sep = "\t",stringsAsFactors = FALSE,header = TRUE)


GO_GSEA_dn <- GO_GSEA_dn[,c(1,6,8)]

GO_GSEA_dn_sig <- GO_GSEA_dn[which(GO_GSEA_dn$FDR.q.val < 0.25 & GO_GSEA_dn$NES < (-1)),]

GO_GSEA_dn_sig$FDR.q.val <- -log10(GO_GSEA_dn_sig_select$FDR.q.val )
 

GO_GSEA_dn_sig_select <- GO_GSEA_dn_sig%>% filter(grepl("IMMUN",NAME)) 

GO_GSEA_dn_sig_select$FDR.q.val <- -log10(GO_GSEA_dn_sig_select$FDR.q.val )



theTable_GO <- within(GO_GSEA_dn_sig, 
                          Pathways <- factor(NAME,levels=GO_GSEA_dn_sig[order(GO_GSEA_dn_sig$FDR.q.val,decreasing=F),]$NAME))



theTable_immune <- within(GO_GSEA_dn_sig_select, 
                   Pathways <- factor(NAME,levels=GO_GSEA_dn_sig_select[order(GO_GSEA_dn_sig_select$FDR.q.val,decreasing=F),]$NAME))


png(filename="final_result/GO_BP_DN_REG.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

ggplot(data=theTable_GO[1:10,], aes(y=NES, x=NAME)) +
  geom_bar(stat="identity",aes(fill = FDR.q.val))+
  scale_fill_gradient(low ="#CC79A7", high = "#D55E00")+
  ylab("Normalized Enrichment Score")+
  xlab("")+
  geom_hline(yintercept = -1.0,linetype='dashed')+
  theme_few()+
  theme(axis.text.x = element_text(face = "bold",angle = 75,vjust = 1,hjust = 1,size =8,family = "arial",colour = "black"),                      # adjusting the position
        axis.title.y = element_text(face = "bold",size = 8,family = "arial"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 8,family = "arial",colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size = 10))


dev.off()

################################



png(filename="final_result/stromalscore2.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(estimate_score, aes(x=Condition, y=StromalScore, fill=Condition)) + 
  geom_boxplot()+
  scale_y_continuous(name = 'StromalScore') +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "t.test")


dev.off()


png(filename="final_result/immunescore2.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(estimate_score, aes(x=Condition, y=ImmuneScore, fill=Condition)) + 
  geom_boxplot()+
  scale_y_continuous(name = 'ImmuneScore') +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "t.test")

dev.off()


png(filename="final_result/estimatescore2.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(estimate_score, aes(x=Condition, y=ESTIMATEScore, fill=Condition)) + 
  geom_boxplot()+
  scale_y_continuous(name = 'ESTIMATEScore') +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "t.test")


dev.off()



png(filename="final_result/TumorPurity2.png",width=800,height=600,unit="px",bg="transparent",res = 150)


ggplot(estimate_score, aes(x=Condition, y=TumorPurity, fill=Condition)) + 
  geom_boxplot()+
  scale_y_continuous(name = 'TumorPurity') +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9),
        axis.line = element_line(size=0.5, colour = "black"))+
  scale_fill_manual(values = c("#2166AC","#C11C0F"))+ 
  stat_compare_means(method = "t.test")



dev.off()
