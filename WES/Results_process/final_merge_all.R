
MLS_maf = read.maf(maf = "MLS.WES.merge.maf")

png(filename="final_result_MUT/MLS_maf_summary.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

plotmafSummary(maf= MLS_maf, rmOutlier = TRUE, addStat = 'median',top = 10,textSize = 0.5,showBarcodes = TRUE)

dev.off()

library(maftools)

png(filename="final_result_MUT/MLS_maf_oncoplot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

oncoplot(maf = MLS_maf, top = 30, showTumorSampleBarcodes = T,fontSize = 0.5,removeNonMutated=FALSE,writeMatrix = TRUE)

dev.off()

RMLS_maf = read.maf(maf = "RMLS.WES.merge.maf")


RMLS_maf

library(openxlsx)

library(data.table)

select_maf <- subset(RMLS_maf_dt,RMLS_maf_dt$Hugo_Symbol== "MICAL3")

MICAL3_alpha_missense <- fread("AlphaMissense-Search-Q7RTP6.tsv")

MICAL3_alpha_missense2 <- MICAL3_alpha_missense[,c(3,5,6,10:13)]

names(MICAL3_alpha_missense2)[1:2] <- c("Hugo_Symbol","aaChange")


aa3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His","Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
aa1 <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K",
         "M", "F", "P", "S", "T", "W", "Y", "V")

library(gsubfn)

select_maf$aaChange <- gsubfn('[A-Z]', setNames(as.list(aa3), aa1), select_maf$aaChange)

library(openxlsx)

MICAL3_alpha_mut <- merge(MICAL3_alpha_missense2,select_maf,"aaChange")

write.xlsx(select_maf,"MICAL3_mutation_RMLS.xlsx")


write.xlsx(MICAL3_alpha_mut,"MICAL3_mutation_site_RMLS.xlsx")


png(filename="final_result_MUT/RMLS_maf_summary.png",width=1000,height=800,unit="px",bg="transparent",res=150)

plotmafSummary(maf= RMLS_maf, rmOutlier = TRUE, addStat = 'median',top = 10,textSize = 0.5,showBarcodes = TRUE)

dev.off()


png(filename="final_result_MUT/RMLS_maf_oncoplot.png",width=1000,height=800,unit="px",bg="transparent",res=150)


oncoplot(maf = RMLS_maf, top = 30, showTumorSampleBarcodes = T,fontSize = 0.5)

dev.off()


#Hippopathway

oncoplot(maf = MLS_maf , pathways = "sigpw", gene_mar = 8, fontSize = 0.6, topPathways = 5)


oncoplot(maf = RMLS_maf , pathways = "sigpw", gene_mar = 8, fontSize = 0.6, topPathways = 5)

MLS.vs.RMLS <- mafCompare(m1 = MLS_maf, m2 = RMLS_maf, m1Name = 'MLS', m2Name = 'RMLS', minMut = 5)



png(filename="final_result_MUT/compare_forest.png",width=800,height=600,unit="px",bg="transparent",res = 150)

forestPlot(mafCompareRes = MLS.vs.RMLS , pVal = 0.1,fdr = 0.25,color = c("#C11C0F","#2166AC"),lineWidth = 1.5,geneFontSize = 1)

dev.off()

###############################################################

lollipopPlot(maf = RMLS_maf, gene = "MICAL3", AACol = 'aaChange', showMutationRate = TRUE)


lollipopPlot(maf = RMLS_maf, gene = "CDH2", AACol = 'aaChange', showMutationRate = TRUE)

lollipopPlot(maf = RMLS_maf, gene = "MDM2", AACol = 'aaChange', showMutationRate = TRUE)

###############################################################################################
RMLS.titv = titv(maf = RMLS_maf, plot = FALSE, useSyn = TRUE)

MLS.titv = titv(maf = MLS_maf, plot = FALSE, useSyn = TRUE)

snvcolors = c( "#C11C0F","#2166AC","#D55E00","#CC79A7","#1b9e77","#e6ab02" )



names(snvcolors) = c("C>A", "C>T", "T>C", "T>A", "C>G", "T>G")

#plot titv summary


png(filename="final_result_MUT/RMLS_titv.png",width=800,height=600,unit="px",bg="transparent",res=100)

plotTiTv(res = RMLS.titv,plotType = 'box',color = snvcolors,baseFontSize = 1,axisTextSize = c(1,1))

dev.off()



png(filename="final_result_MUT/MLS_titv.png",width=800,height=600,unit="px",bg="transparent",res=100)

plotTiTv(res = MLS.titv,plotType = 'box',color = snvcolors,baseFontSize = 1,axisTextSize = c(1,1))

dev.off()


plotTiTv(res = MLS.titv,plotType = 'box',color = snvcolors)

plotmafSummary(maf = MLS_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#######################################################################################################
library(BSgenome)

library(maftools)


MLS_maf.tnm = trinucleotideMatrix(maf = MLS_maf,
                                      ref_genome = "BSgenome.Hsapiens.UCSC.hg38")


library(NMF)



MLS_maf.sig.ext <- extractSignatures(mat = MLS_maf.tnm, 
                                         n = 3,
                                         pConstant = 0.1,
                                         parallel = 1)




png(filename="final_result_MUT/MLS_signature_sig.png",width=1000,height=800,unit="px",bg="transparent",res=150)

plotSignatures(nmfRes = MLS_maf.sig.ext , 
               title_size = 1.2,
               contributions = FALSE,
               show_title = TRUE,
               sig_db = 'legacy')

dev.off()


#################################################################################################

RMLS_maf.tnm = trinucleotideMatrix(maf = RMLS_maf,
                                  ref_genome = "BSgenome.Hsapiens.UCSC.hg38")


library(NMF)



RMLS_maf.sig.ext <- extractSignatures(mat = RMLS_maf.tnm, 
                                     n = 3,
                                     pConstant = 0.1,
                                     parallel = 1)


plotSignatures(RMLS_maf.sig.ext)

png(filename="final_result_MUT/RMLS_signature_sig.png",width=1000,height=800,unit="px",bg="transparent",res=150)

plotSignatures(nmfRes = RMLS_maf.sig.ext, 
               title_size = 1.2,
               contributions = FALSE,
               show_title = TRUE,
               sig_db = 'legacy')

dev.off()
########################################################################


MLS.sig = extractSignatures(mat = MLS_maf.tnm , n =3)

MLS.sig.cosm = compareSignatures(nmfRes = MLS.sig , sig_db = "legacy")

MLS.sig.cosm_SBS = compareSignatures(nmfRes = MLS.sig, sig_db = "SBS")

RMLS.sig = extractSignatures(mat = RMLS_maf.tnm , n = 3)

RMLS.sig.cosm = compareSignatures(nmfRes = RMLS.sig , sig_db = "legacy")

RMLS.sig.cosm_SBS = compareSignatures(nmfRes = RMLS.sig, sig_db = "SBS")


library('pheatmap')
pheatmap::pheatmap(mat = MLS.sig.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

pheatmap::pheatmap(mat = RMLS.sig.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

#########################################################################
library(dndscv)

MLS_maf_dt <- MLS_maf@data

MLS_mut <- MLS_maf_dt[,c(6,1,2,4,5)]

colnames(MLS_mut) <- c("sampleID","chr","pos","ref","mut")

MLS_mut$chr = gsub("chr","",as.vector(MLS_mut$chr))

MLS_mut <- na.omit(MLS_mut)

MLS_mut <- unique(MLS_mut)

MLS_dndscv <- dndscv(MLS_mut,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,refdb = "RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv = NULL)

MLS_sel_cv = MLS_dndscv$sel_cv

signif_genes_MLS = MLS_sel_cv[MLS_sel_cv$qglobal_cv<0.05, c("gene_name","qglobal_cv")]

RMLS_maf_dt <- RMLS_maf@data

RMLS_mut <- RMLS_maf_dt[,c(6,1,2,4,5)]

colnames(RMLS_mut) <- c("sampleID","chr","pos","ref","mut")

RMLS_mut$chr = gsub("chr","",as.vector(RMLS_mut$chr))

RMLS_dndscv <- dndscv(RMLS_mut,max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,refdb = "RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv = NULL)

RMLS_sel_cv = RMLS_dndscv$sel_cv

RMLS_signif_genes = RMLS_sel_cv[RMLS_sel_cv$qglobal_cv<0.05, c("gene_name","qglobal_cv")]

library(dplyr)

RMLS_signif_only_genes <- anti_join(RMLS_signif_genes ,signif_genes_MLS,"gene_name")

library(openxlsx)

write.xlsx(RMLS_signif_only_genes,"RMLS_significant_only_genes.xlsx")

MLS_signif_only_genes <- anti_join(signif_genes_MLS,RMLS_signif_genes ,"gene_name")
################################################################################


mls.mutsig.corrected = prepareMutSig(maf = MLS_maf,fn = "MLS" )



rmls.mutsig.corrected = prepareMutSig(maf = RMLS_maf,fn = "RMLS" )

#./run_MutSigCV.sh ~/MatlabMCR/v901/ ../MLS.mutSig.maf ../ref/exome_full192.coverage.txt ../ref/gene.covariates.txt MLS_mutsig ../ref/mutation_type_dictionary_file.txt ../ref/chr_files_hg38.txt


###########################################################################################
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


Enriched_go <- enrichr(genes = RMLS_signif_only_genes$gene_name, databases = dbs_go)


Enriched_go_mls <- enrichr(genes = MLS_signif_only_genes$gene_name, databases = dbs_go)

png(filename="final_result_MUT/RMLS_sig_genes_go_bp.png",width=1000,height=800,unit="px",bg="transparent",res=100)

plotEnrich(Enriched_go[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for driver mutation in RMLS")

dev.off()

plotEnrich(Enriched_go_mls[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

Enriched_pw <- enrichr(genes = RMLS_signif_only_genes$gene_name, databases = dbs_pw)

Enriched_pw_mls <- enrichr(genes = MLS_signif_only_genes$gene_name, databases = dbs_pw)

plotEnrich(Enriched_pw[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")


png(filename="final_result_MUT/RMLS_sig_genes_kegg.png",width=1000,height=800,unit="px",bg="transparent",res=100)

plotEnrich(Enriched_pw[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for driver mutation in RMLS")

dev.off()

plotEnrich(Enriched_pw_mls[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

plotEnrich(Enriched_pw_mls[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value")

dbs_hm <-  c("MSigDB_Hallmark_2020","MSigDB_Oncogenic_Signatures")

Enriched_hm <- enrichr(genes =  RMLS_signif_only_genes$gene_name, databases =dbs_hm)


Enriched_hm_mls <- enrichr(genes =  MLS_signif_only_genes$gene_name, databases =dbs_hm)

plotEnrich(Enriched_hm[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

plotEnrich(Enriched_hm_mls[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")


Enriched_ppi <- enrichr(genes = RMLS_signif_only_genes$gene_name, databases = dbs_ppi)


plotEnrich(Enriched_ppi[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "PPI enrichment analysis for driver mutation in RMLS")


Enriched_ppi_mls <- enrichr(genes = MLS_signif_only_genes$gene_name, databases = dbs_ppi)

ppi_MLS <- Enriched_ppi_mls$PPI_Hub_Proteins



ppi_RMLS <- Enriched_ppi$PPI_Hub_Proteins

plotEnrich(Enriched_ppi_mls[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")


only_ppi_RMLS <- anti_join(ppi_RMLS,ppi_MLS,"Term")


png(filename="final_result_MUT/MLS_somatic_interaction.png",width=1000,height=800,unit="px",bg="transparent",res=150)

somaticInteractions(maf = MLS_maf, top = 30, pvalue = c(0.05, 0.1),fontSize = 0.5,)

dev.off()


png(filename="final_result_MUT/RMLS_somatic_interaction.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


somaticInteractions(maf = RMLS_maf, top = 30, pvalue = c(0.05, 0.1),fontSize = 0.5)

dev.off()
####################################################################

write.mafSummary(MLS_1_maf,basename = "MLS_1")

MLS_1_maf <- subsetMaf(MLS_maf,tsb = "MLS_1")


write.mafSummary(RMLS_1_maf,basename = "RMLS_1")

RMLS_1_maf <- subsetMaf(RMLS_maf,tsb = "RMLS_1")


####################################################################

MLS_all_maf <- read.maf(maf = "MLS.WES.merge.maf",gisticAllLesionsFile =  'CNVkit_gistic_MLS/all_lesions.conf_99.txt',
                         gisticAmpGenesFile = 'CNVkit_gistic_MLS/amp_genes.conf_99.txt',
                         gisticDelGenesFile = 'CNVkit_gistic_MLS/del_genes.conf_99.txt',
                         gisticScoresFile = 'CNVkit_gistic_MLS/scores.gistic')


library(maftools)

library(MesKit)

MLS.gistic = readGistic(gisticAllLesionsFile =  'CNVkit_gistic_MLS/all_lesions.conf_99.txt',
                        gisticAmpGenesFile = 'CNVkit_gistic_MLS/amp_genes.conf_99.txt',
                        gisticDelGenesFile = 'CNVkit_gistic_MLS/del_genes.conf_99.txt',
                        gisticScoresFile = 'CNVkit_gistic_MLS/scores.gistic')


png(filename="final_result_MUT/MLS_gistic_chrom_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

gisticChromPlot(gistic = MLS.gistic , ref.build = "hg38",markBands = "all",mutGenes = TRUE,fdrCutOff = 0.25,txtSize = 0.8,cytobandTxtSize = 0.5,y_lims = c(-1.0,1.0))

dev.off()


png(filename="final_result_MUT/MLS_gistic_Onco_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

gisticOncoPlot(MLS.gistic,sortByAnnotation = TRUE,showTumorSampleBarcodes = TRUE,annotationFontSize = 5,colors = c(Amp = "red", Del = "blue"),fontSize =0.5,top = 30)

dev.off()


library(tidyverse)

(MLS_sig_cytoband <- MLS.gistic@cytoband.summary %>% filter(qvalues<0.25) %>% .$Unique_Name)


MLS_sig_AP_cytoband <- MLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  .$Unique_Name %>% 
  .[grepl(pattern = "AP",.)]

MLS_sig_AP_gene <- MLS.gistic@data %>% filter(Cytoband %in% MLS_sig_AP_cytoband) %>% .$Hugo_Symbol 

MLS_sig_AP_gene <- unique(MLS_sig_AP_gene)

MLS_sig_DP_cytoband <- MLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  .$Unique_Name %>% 
  .[grepl(pattern = "DP",.)]

MLS_sig_DP_gene <- MLS.gistic@data %>% filter(Cytoband %in% MLS_sig_DP_cytoband) %>% .$Hugo_Symbol 

MLS_sig_DP_gene <- unique(MLS_sig_DP_gene)

library(ggplot2)

library(enrichR)

dbs <- listEnrichrDbs()

dbs <- dbs[order(dbs$libraryName),]

dbs$libraryName

dbs_go <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")

dbs_go

dbs_pw <- c("KEGG_2021_Human", "WikiPathways_2024_Human", "Reactome_2022")

dbs_dd <- c("PheWeb_2019", "ClinVar_2019")

dbs_ppi <- c("PPI_Hub_Proteins")


dbs_hm <-  c("MSigDB_Hallmark_2020","MSigDB_Oncogenic_Signatures")

Enriched_MLS_AP_go <- enrichr(genes = MLS_sig_AP_gene , databases = dbs_go)


Enriched_MLS_AP_pw <- enrichr(genes = MLS_sig_AP_gene , databases = dbs_pw)

Enriched_MLS_DP_go <- enrichr(genes = MLS_sig_DP_gene , databases = dbs_go)


Enriched_MLS_DP_pw <- enrichr(genes = MLS_sig_DP_gene , databases = dbs_pw)


plotEnrich(Enriched_MLS_DP_pw[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for genes of deletion in MLS")

png(filename="FINAL_MUT/final_result_MUT/MLS_AP_GO_BP_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_MLS_AP_go[[3]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for genes of amplification in MLS")

dev.off()


png(filename="FINAL_MUT/final_result_MUT/MLS_DP_GO_BP_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_MLS_DP_go[[3]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for genes of deletion in MLS")

dev.off()
#################################################################################
RMLS_all_maf <- read.maf(maf = "RMLS.WES.merge.maf",gisticAllLesionsFile =  'CNVkit_gistic_RMLS/all_lesions.conf_99.txt',
                         gisticAmpGenesFile = 'CNVkit_gistic_RMLS/amp_genes.conf_99.txt',
                         gisticDelGenesFile = 'CNVkit_gistic_RMLS/del_genes.conf_99.txt',
                         gisticScoresFile = 'CNVkit_gistic_RMLS/scores.gistic')

RMLS.gistic = readGistic(gisticAllLesionsFile =  'CNVkit_gistic_RMLS/all_lesions.conf_99.txt',
                        gisticAmpGenesFile = 'CNVkit_gistic_RMLS/amp_genes.conf_99.txt',
                        gisticDelGenesFile = 'CNVkit_gistic_RMLS/del_genes.conf_99.txt',
                        gisticScoresFile = 'CNVkit_gistic_RMLS/scores.gistic')


png(filename="final_result_MUT/RMLS_gistic_chrom_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

gisticChromPlot(gistic = RMLS.gistic , ref.build = "hg19",markBands = c("1p36.32","6q26","19q13.43"),fdrCutOff = 0.5,txtSize = 0.5,cytobandTxtSize = 0.5,y_lims = c(-1.0,1.0),mutGenes = c("HES5" ,
"THBS2","ZNF460") ,mutGenesTxtSize = 1.0,maf = RMLS_all_maf,cytobandOffset = 0.01 )

dev.off()


png(filename="final_result_MUT/RMLS_gistic_Onco_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

gisticOncoPlot(RMLS.gistic,sortByAnnotation = TRUE,showTumorSampleBarcodes = TRUE,annotationFontSize = 5,colors = c(Amp = "red", Del = "blue"),fontSize =0.5,top = 30)

dev.off()


(RMLS_sig_cytoband <- RMLS.gistic@cytoband.summary %>% filter(qvalues<0.25) %>% .$Unique_Name)


RMLS_sig_AP_cytoband <- RMLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  .$Unique_Name %>% 
  .[grepl(pattern = "AP",.)]

RMLS_sig_AP_gene <- RMLS.gistic@data %>% filter(Cytoband %in% RMLS_sig_AP_cytoband) %>% .$Hugo_Symbol 

RMLS_sig_AP_gene <- unique(RMLS_sig_AP_gene)

RMLS_sig_DP_cytoband <- RMLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  .$Unique_Name %>% 
  .[grepl(pattern = "DP",.)]

RMLS_sig_DP_cytoband 

MLS_sig_DP_cytoband 

RMLS.gistic@numericMatrix

RMLS_sig_DP_gene <- RMLS.gistic@data %>% filter(Cytoband %in% RMLS_sig_DP_cytoband) %>% .$Hugo_Symbol 

RMLS_sig_DP_gene <- unique(RMLS_sig_DP_gene)


Enriched_RMLS_AP_go <- enrichr(genes = RMLS_sig_AP_gene , databases = dbs_go)

Enriched_RMLS_DP_go <- enrichr(genes = RMLS_sig_DP_gene , databases = dbs_go)


Enriched_RMLS_AP_pw <- enrichr(genes = RMLS_sig_AP_gene , databases = dbs_pw)

Enriched_RMLS_DP_pw <- enrichr(genes = RMLS_sig_DP_gene , databases = dbs_pw)

png(filename="FINAL_MUT/final_result_MUT/RMLS_AP_GO_BP_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_RMLS_AP_go[[3]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for genes of amplification in RMLS")

dev.off()


plotEnrich(Enriched_RMLS_AP_pw[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for genes of amplification in RMLS")


png(filename="final_result_MUT/RMLS_DP_GO_BP_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_RMLS_DP_go[[3]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for genes of deletion in RMLS")

dev.off()

Enriched_RMLS_DP_go[[3]]

png(filename="final_result_MUT/RMLS_DP_KEGG_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_RMLS_DP_pw[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for genes of deletion in RMLS")

dev.off()

RMLS_sig_DP_gene_df = data.frame(Gene = RMLS_sig_DP_gene,RMLS ="Del")


RMLS_sig_AP_gene_df = data.frame(Gene = RMLS_sig_AP_gene,RMLS = "Amp")


MLS_sig_DP_gene_df = data.frame(Gene = MLS_sig_DP_gene,MLS = "Del")

MLS_sig_AP_gene_df = data.frame(Gene = MLS_sig_AP_gene, MLS = "Amp")


RMLS_sig_cnv_genes <- rbind(RMLS_sig_AP_gene_df,RMLS_sig_DP_gene_df)


MLS_sig_cnv_genes <- rbind(MLS_sig_AP_gene_df,MLS_sig_DP_gene_df)

all_cnv_sig <- merge(RMLS_sig_cnv_genes,MLS_sig_cnv_genes)


library(dplyr)

RMLS_only_DP_gene <- anti_join(RMLS_sig_DP_gene_df,MLS_sig_DP_gene_df)


Enriched_only_RMLS_DP_pw <- enrichr(genes = RMLS_only_DP_gene$Gene , databases = dbs_pw)


Enriched_only_RMLS_DP_go <- enrichr(genes = RMLS_only_DP_gene$Gene , databases = dbs_go)


png(filename="final_result_MUT/ONLY_RMLS_DP_GO_BP_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_only_RMLS_DP_go[[3]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for genes of deletion in RMLS")

dev.off()

png(filename="final_result_MUT/ONLY_RMLS_DP_KEGG_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_only_RMLS_DP_pw[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for genes of deletion in RMLS")

dev.off()

###################################################################
MLS.sig = oncodrive(maf = MLS_maf, AACol = 'aaChange', minMut = 1, pvalMethod = 'zscore')

MLS.select.sig <- MLS.sig[,c(1,19)]

names(MLS.select.sig) <- c("Hugo_Symbol","Oncodrive")

MLS.select.sig<- subset(MLS.select.sig,MLS.select.sig$Oncodrive < 0.05)

signif_genes_MLS_select <- signif_genes_MLS

names(signif_genes_MLS_select) <- c("Hugo_Symbol","dndscv")

MLS_onco_dnd_sig <- merge(MLS.select.sig,signif_genes_MLS_select)


RMLS.sig = oncodrive(maf = RMLS_maf, AACol = 'aaChange', minMut = 1, pvalMethod = 'zscore')


RMLS.select.sig <- RMLS.sig[,c(1,19)]

names(RMLS.select.sig) <- c("Hugo_Symbol","Oncodrive")

RMLS.select.sig<- subset(RMLS.select.sig,RMLS.select.sig$Oncodrive < 0.05)


signif_genes_RMLS_select <- RMLS_signif_genes

names(signif_genes_RMLS_select) <- c("Hugo_Symbol","dndscv")

RMLS_onco_dnd_sig <- merge(RMLS.select.sig,signif_genes_RMLS_select)

library(openxlsx)

write.xlsx(MLS_onco_dnd_sig,"MLS_onco_dndscv_sig.xlsx")


write.xlsx(RMLS_onco_dnd_sig,"RMLS_onco_dndscv_sig.xlsx")


plotVaf(maf = RMLS_maf, vafCol = 'VAF',genes = RMLS_onco_dnd_sig$"Hugo_Symbol")


plotVaf(maf = MLS_maf, vafCol = 'VAF',genes = MLS_onco_dnd_sig$"Hugo_Symbol")

library(dplyr)


MLS_mutation <- subset(MLS_maf@data,MLS_maf@data$Hugo_Symbol %in% MLS_onco_dnd_sig$Hugo_Symbol)

MLS_mutation <- MLS_mutation[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification","aaChange","VAF" )]


MLS_mutation <- merge(MLS_mutation,MLS_onco_dnd_sig)


RMLS_mutation <- subset(RMLS_maf@data,RMLS_maf@data$Hugo_Symbol %in% RMLS_onco_dnd_sig$Hugo_Symbol)

RMLS_mutation <- RMLS_mutation[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification","aaChange","VAF" )]

RMLS_mutation <- merge(RMLS_mutation,RMLS_onco_dnd_sig)

all_sig_mut_genes <- rbind(data.frame(Hugo_Symbol = unique(MLS_mutation$Hugo_Symbol)),data.frame(Hugo_Symbol =unique(RMLS_mutation$Hugo_Symbol)))

oncoplot(maf = RMLS_maf,genes =c(all_sig_mut_genes$Hugo_Symbol),keepGeneOrder = TRUE)


oncoplot(maf = MLS_maf,genes =c(all_sig_mut_genes$Hugo_Symbol),keepGeneOrder = TRUE)


summary_mls_genes <- MLS_maf@gene.summary

top10_mls_genes <-  summary_mls_genes$Hugo_Symbol[1:12]

top10_mls_genes 

summary_rmls_genes <- RMLS_maf@gene.summary

top10_rmls_genes <-  summary_rmls_genes$Hugo_Symbol[1:12]

top10_rmls_genes



oncoplot(maf = RMLS_maf,genes =c(top10_mls_genes, top10_rmls_genes) ,keepGeneOrder = TRUE)


oncoplot(maf = RMLS_all_maf,genes =c("HES5","THBS2","ZNF460") ,keepGeneOrder = TRUE)


oncoplot(maf = MLS_all_maf,genes =c("HES5","THBS2","ZNF460") ,keepGeneOrder = TRUE)



oncoplot(maf = RMLS_all_maf,top = 30 ,keepGeneOrder = TRUE,MAT)

oncoplot(maf = MLS_maf,genes =c(top10_mls_genes, top10_rmls_genes) ,keepGeneOrder = TRUE)


library(ComplexHeatmap)

MLS_mut_mat <- read.table("MLS_top_onco_matrix.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)

MLS_mut_mat$MLS_5 <- gsub(x = MLS_mut_mat$MLS_5 ,pattern = 0,replacement = "")


RMLS_mut_mat <- read.table("RMLS_top_onco_matrix.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)

  


top10_each_genes_vaf_MLS = subsetMaf(maf =MLS_maf, genes = c(top10_mls_genes, top10_rmls_genes), fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]


top10_each_genes_vaf_RMLS = subsetMaf(maf =RMLS_maf, genes = c(top10_mls_genes, top10_rmls_genes), fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]


#MLS_mut_vaf <- MLS_mutation[,mean(VAF, na.rm = TRUE), Hugo_Symbol]

library(tibble)

#MLS_mut_vaf <- column_to_rownames(MLS_mut_vaf,"Hugo_Symbol")

names(MLS_mut_vaf) = "VAF"

library(ComplexHeatmap)

MLS_mut_mat_df <- rownames_to_column(MLS_mut_mat,"Hugo_Symbol")
 
top20_mut_df <- data.frame(Hugo_Symbol=MLS_mut_mat_df$Hugo_Symbol)

top20_MLS_vaf <- merge(top20_mut_df,top10_each_genes_vaf_MLS,sort = FALSE)

colnames(top20_MLS_vaf)[2] <- "VAF"

top20_RMLS_vaf <- merge(top20_mut_df,top10_each_genes_vaf_RMLS,sort = FALSE,all.x = TRUE)

colnames(top20_RMLS_vaf)[2] <- "VAF"

top20_RMLS_vaf[is.na(top20_RMLS_vaf)] = 0


top20_RMLS_vaf <- merge(top20_mut_df,top20_RMLS_vaf,sort = FALSE,all.x = TRUE)
############################################################

RMLS_sig_res_genes <-  RMLS_sel_cv[,c("gene_name","qglobal_cv")]

names(RMLS_sig_res_genes) <- c("Hugo_Symbol","q_value")

RMLS_sig_res_genes$q_value <- -log10(RMLS_sig_res_genes$q_value)


top20_RMLS_vaf <- merge(top20_RMLS_vaf,RMLS_sig_res_genes,sort = FALSE,all.x = TRUE)

MLS_sig_res_genes <- MLS_sel_cv[,c("gene_name","qglobal_cv")]


names(MLS_sig_res_genes) <- c("Hugo_Symbol","q_value")

MLS_sig_res_genes$q_value <- -log10(MLS_sig_res_genes$q_value)

top20_MLS_vaf <- merge(top20_MLS_vaf,MLS_sig_res_genes,sort = FALSE)

###############################################################################

library(ComplexHeatmap)

library(ggthemes)
library(ggplot2)

col <- c("Nonsense_Mutation" = "red","Missense_Mutation" = "#008000","In_Frame_Del"="yellow","Frame_Shift_Del"="blue","Translation_Start_Site" ="orange","Multi_Hit"="black")


col2 <- c("Nonsense_Mutation" = "red","Missense_Mutation" = "#008000","In_Frame_Del"="yellow","Frame_Shift_Del"="blue","Multi_Hit"="black")



alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  }
)



alter_fun2 <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  }
)


heatmap_legend_param <- list(title = "", at = c("Nonsense_Mutation" , "Missense_Mutation","In_Frame_Del","Frame_Shift_Del","Translation_Start_Site","Multi_Hit"),labels = c("Nonsense_Mutation","Missense_Mutation", 
                                                                                                                                                                            "In_Frame_Del","Frame_Shift_Del","Translation_Start_Site","Multi_Hit"),labels_gp = gpar(fontsize = 5.8,fontface = 2,fontfamily="Arial"),ncol = 3, by_row = TRUE)



heatmap_legend_param2 <- list(title = "", at = c("Nonsense_Mutation" , "Missense_Mutation","In_Frame_Del",
                                                 "Frame_Shift_Del","Multi_Hit"), 
                              labels = c("Nonsense_Mutation","Missense_Mutation", 
                                         "In_Frame_Del","Frame_Shift_Del","Multi_Hit"))

h1 = rowAnnotation("q-value" = anno_barplot(top20_MLS_vaf$q_value,axis_param=list(direction = "reverse",gp=gpar(fontsize = 6,fontface = 2,fontfamily="Arial"),at = c(0, 1, 2))),width = unit(1.0, "cm"), annotation_name_gp = gpar(fontsize=6,fontface = 2,fontfamily="Arial")) + oncoPrint(MLS_mut_mat ,show_pct = FALSE,right_annotation =rowAnnotation(Type= c("Both",rep("MLS",5),"Both",rep("MLS",5),rep("RMLS",10)), col = list(Type = c("MLS" = "#2166AC","RMLS" = "#C11C0F" ,"Both" = "#1B9E77")),annotation_legend_param = list(Type = list(title = "",labels_gp = gpar(fontsize = 5.8,fontface = 2,fontfamily="Arial"),ncol = 3, by_row = TRUE)), width = unit(1, "mm"), annotation_name_gp = gpar(fontsize=6,fontface = 2,fontfamily="Arial")),left_annotation = rowAnnotation("VAF" = anno_barplot(top20_MLS_vaf$VAF,axis_param = list(direction = "reverse",gp=gpar(fontsize = 6,fontface = 2,fontfamily="Arial"),at = c(0, 0.2, 0.4))),width = unit(1.0, "cm"),annotation_name_gp = gpar(fontsize=6,fontface = 2,fontfamily="Arial")),
               alter_fun = alter_fun, col = col,
               column_labels = colnames(MLS_mut_mat),
               show_column_names = TRUE,
               remove_empty_columns = FALSE,
               remove_empty_rows = FALSE,top_annotation = NULL,
               heatmap_legend_param = heatmap_legend_param,
               column_names_gp=gpar(fontsize = 6,fontface = 2,fontfamily="Arial"),
               row_names_gp=gpar(fontsize =6,fontface = 2,fontfamily="Arial"))  

h2 = oncoPrint(RMLS_mut_mat,show_row_names = FALSE, show_pct = FALSE, right_annotation = rowAnnotation("VAF" = anno_barplot(top20_RMLS_vaf$VAF,axis_param=list(gp=gpar(fontsize = 6,fontface = 2,fontfamily="Arial"),at = c(0, 0.2, 0.4)),width = unit(1.0, "cm")),                                                             annotation_name_gp = gpar(fontsize=6,fontface = 2,fontfamily="Arial")),
                alter_fun = alter_fun2, col = col2, column_labels = colnames(RMLS_mut_mat) ,show_column_names = TRUE,
               left_annotation = NULL,
                remove_empty_columns = FALSE, top_annotation = NULL,
                remove_empty_rows = FALSE ,row_labels = NULL,show_heatmap_legend = FALSE,column_names_gp=gpar(fontsize = 6,fontface = 2,fontfamily="Arial")) + rowAnnotation("q-value" = anno_barplot(top20_RMLS_vaf$q_value,axis_param=list(gp=gpar(fontsize = 6,fontface = 2,fontfamily="Arial"),at = c(0, 1, 2))),width = unit(1.0, "cm"), annotation_name_gp = gpar(fontsize=6,fontface = 2,fontfamily="Arial"))


png(filename="final_result_MUT/coonco_plot.png",width=800,height=600,unit="px",bg="transparent",res = 150)

draw(h1+ h2,heatmap_legend_side = "bottom")

dev.off()
##############################################################################


RMLS_sig_mutation_vaf = subsetMaf(maf =RMLS_maf, genes =c(all_sig_mut_genes$Hugo_Symbol), fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]


MLS_sig_mutation_vaf = subsetMaf(maf =MLS_maf, genes =c(all_sig_mut_genes$Hugo_Symbol), fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]

RMLS_mutation_re <- RMLS_mutation[,c("Hugo_Symbol","VAF","dndscv","Oncodrive")]

RMLS_mutation_re <- unique(RMLS_mutation_re)



MLS_sig_mut_mat <- read.table("MLS_sig_onco_matrix.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)

MLS_sig_mut_mat$MLS_5 <- gsub(x = MLS_sig_mut_mat$MLS_5 ,pattern = 0,replacement = "")

MLS_sig_mut_mat_df <- rownames_to_column(MLS_sig_mut_mat,"Hugo_Symbol")

all_sig_mut_df <- data.frame(Hugo_Symbol=MLS_sig_mut_mat_df$Hugo_Symbol)

all_sig_MLS_vaf <- merge(all_sig_mut_df,MLS_sig_mutation_vaf ,sort = FALSE,all.x = TRUE)

colnames(all_sig_MLS_vaf)[2] <- "VAF"

all_sig_MLS_vaf[is.na(all_sig_MLS_vaf)] = 0


all_sig_MLS_vaf <- merge(all_sig_mut_df,all_sig_MLS_vaf,sort = FALSE,all.x = TRUE)

all_dndscv_genes_MLS = MLS_sel_cv[, c("gene_name","qglobal_cv")]

names(all_dndscv_genes_MLS) <- c("Hugo_Symbol","dndscv")

all_sig_MLS_vaf2 <- merge(all_sig_MLS_vaf,all_dndscv_genes_MLS,sort = FALSE,all.x = TRUE)


all_MLS.select.sig <- MLS.sig[,c(1,19)]

names(all_MLS.select.sig) <- c("Hugo_Symbol","Oncodrive")

all_sig_MLS_vaf2 <- merge(all_sig_MLS_vaf2,all_MLS.select.sig,sort = FALSE,all.x = TRUE)


all_sig_MLS_vaf2 [is.na(all_sig_MLS_vaf2 )] = 0


all_sig_MLS_vaf2$dndscv <- -log10(all_sig_MLS_vaf2$dndscv)

all_sig_MLS_vaf2
###########################################################################

RMLS_sig_mut_mat <- read.table("RMLS_sig_onco_matrix.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)


RMLS_sig_mut_mat_df <- rownames_to_column(RMLS_sig_mut_mat,"Hugo_Symbol")

all_sig_RMLS_vaf <- merge(all_sig_mut_df,RMLS_sig_mutation_vaf ,sort = FALSE,all.x = TRUE)

colnames(all_sig_RMLS_vaf)[2] <- "VAF"

all_sig_RMLS_vaf[is.na(all_sig_RMLS_vaf)] = 0


all_sig_RMLS_vaf <- merge(all_sig_mut_df,all_sig_RMLS_vaf,sort = FALSE,all.x = TRUE)


all_dndscv_genes_RMLS = RMLS_sel_cv[, c("gene_name","qglobal_cv")]

names(all_dndscv_genes_RMLS) <- c("Hugo_Symbol","dndscv")

all_sig_RMLS_vaf2 <- merge(all_sig_RMLS_vaf,all_dndscv_genes_RMLS,sort = FALSE,all.x = TRUE)


all_RMLS.select.sig <- RMLS.sig[,c(1,19)]

names(all_RMLS.select.sig) <- c("Hugo_Symbol","Oncodrive")

all_sig_RMLS_vaf2 <- merge(all_sig_RMLS_vaf2,all_RMLS.select.sig,sort = FALSE,all.x = TRUE)


all_sig_RMLS_vaf2 [is.na(all_sig_RMLS_vaf2 )] = 0


all_sig_RMLS_vaf2 <- merge(all_sig_mut_df,all_sig_RMLS_vaf2,sort = FALSE,all.x = TRUE)

all_sig_RMLS_vaf2 

#MLS_mutation_sig <- MLS_mutation[,c("Hugo_Symbol")]


col_sig <- c("Nonsense_Mutation" = "red","Missense_Mutation" = "#008000","In_Frame_Del"="yellow","Frame_Shift_Del"="blue","Nonstop_Mutation" ="skyblue","Multi_Hit"="black")


alter_fun_sig <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col[" Nonstop_Mutation"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  }
)


heatmap_legend_param_sig <- list(title = "", at = c("Nonsense_Mutation" , "Missense_Mutation","In_Frame_Del","Frame_Shift_Del","Nonstop_Mutation","Multi_Hit"),labels = c("Nonsense_Mutation","Missense_Mutation",  "In_Frame_Del","Frame_Shift_Del","Nonstop_Mutation","Multi_Hit"),labels_gp = gpar(fontsize = 5),ncol = 3, by_row = TRUE)

library(RColorBrewer)

library(circlize)

qvalue_col_fun = colorRamp2(c(0, 1, 2), c("white", "grey", "red"))


all_sig_MLS_vaf2$dndscv



h3 =  oncoPrint(MLS_sig_mut_mat ,show_pct = FALSE,right_annotation =rowAnnotation(Type= c(rep("MLS",16),rep("RMLS",17)), col = list(Type = c("MLS" = "blue","RMLS" = "red" )),annotation_legend_param = list(Type = list(title = "",labels_gp = gpar(fontsize = 5),ncol = 3, by_row = TRUE)), width = unit(1, "mm"), annotation_name_gp = gpar(fontsize=5)),left_annotation = rowAnnotation("VAF" = anno_barplot(all_sig_MLS_vaf$VAF,axis_param = list(direction = "reverse",gp=gpar(fontsize = 4),at = c(0, 0.2, 0.4))), width = unit(1.0, "cm"),
 annotation_name_gp = gpar(fontsize=5)),alter_fun = alter_fun_sig, col = col_sig, column_labels = colnames(MLS_sig_mut_mat),show_column_names = TRUE,
             remove_empty_columns = FALSE,remove_empty_rows = FALSE,top_annotation = NULL, heatmap_legend_param = heatmap_legend_param_sig, column_names_gp=gpar(fontsize = 5),row_names_gp=gpar(fontsize = 5))  

h4 = oncoPrint(RMLS_sig_mut_mat,show_row_names = FALSE, show_pct = FALSE, right_annotation = rowAnnotation("VAF" = anno_barplot(all_sig_RMLS_vaf$VAF,axis_param=list(gp=gpar(fontsize = 4),at = c(0, 0.2, 0.4)),width = unit(1.0, "cm")),                                                             annotation_name_gp = gpar(fontsize=5)),
               alter_fun = alter_fun_sig, col = col_sig, column_labels = colnames(RMLS_sig_mut_mat) ,show_column_names = TRUE,
               left_annotation = NULL,
               remove_empty_columns = FALSE, top_annotation = NULL,
               remove_empty_rows = FALSE ,row_labels = NULL,show_heatmap_legend = FALSE,column_names_gp=gpar(fontsize = 5)) 


png(filename="final_result_MUT/coonco_sig_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


draw(h3+ h4,heatmap_legend_side = "bottom")

dev.off()

###############################

RMLS_sig_DP_cytoband_df <- data.frame(DP_cytoband = RMLS_sig_DP_cytoband)

MLS_sig_DP_cytoband_df <- data.frame(DP_cytoband = MLS_sig_DP_cytoband)

library(dplyr)

RMLS_sig_DP_cytoband_df <- strsplit(RMLS_sig_DP_cytoband_df$DP_cytoband, split = ":")

RMLS_sig_DP_cytoband_df <- data.frame(DP_cytoband = unlist(RMLS_sig_DP_cytoband_df))

RMLS_sig_DP_cytoband_df <-  data.frame(DP_cytoband = RMLS_sig_DP_cytoband_df[c(FALSE, TRUE), ])


MLS_sig_DP_cytoband_df <- strsplit(MLS_sig_DP_cytoband_df$DP_cytoband, split = ":")

MLS_sig_DP_cytoband_df <- data.frame(DP_cytoband = unlist(MLS_sig_DP_cytoband_df))

MLS_sig_DP_cytoband_df <-  data.frame(DP_cytoband = MLS_sig_DP_cytoband_df[c(FALSE, TRUE), ])


library(stringr)

rm(list = c("RMLS_all_maf_df"))

RMLS.gistic@cnMatrix

library(dplyr)



RMLS_sig_AP_cytoband_sort <- RMLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband) %>%
  .[grepl( Unique_Name,pattern = "AP",.)]


MLS_sig_AP_cytoband_sort <- MLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband) %>%
  .[grepl( Unique_Name,pattern = "AP",.)]

RMLS_only_sig_AP_cytoband <- anti_join(RMLS_sig_AP_cytoband_sort,MLS_sig_AP_cytoband_sort,"Cytoband")

RMLS_only_sig_AP_cytoband_genes <- RMLS.gistic@data %>% filter(Cytoband %in% RMLS_only_sig_AP_cytoband$Unique_Name) %>% select(Hugo_Symbol, Cytoband) 


RMLS_only_sig_AP_mutation_genes <- merge(RMLS.only.drivermut.sig,RMLS_only_sig_AP_cytoband_genes,"Hugo_Symbol")


RMLS_only_sig_DP_mutation_genes <- merge(RMLS.only.drivermut.sig,RMLS_only_sig_DP_cytoband_genes,"Hugo_Symbol")


RMLS_only_sig_DP_mutation_0.1_genes <- merge(RMLS.only.drivermut.0.1.sig,RMLS_only_sig_DP_cytoband_genes,"Hugo_Symbol")


RMLS_only_sig_DP_mutation_genes <- unique(RMLS_only_sig_DP_mutation_genes)

RMLS_only_sig_DP_mutation_genes

RMLS.only.drivermut.0.1.sig



plotEnrich(Enriched_pw_RMLS_only_DP_genes[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis ")


Enriched_RMLS.only.drivermut.0.1.sig <- enrichr(genes = RMLS.only.drivermut.0.1.sig$Hugo_Symbol, databases = dbs_pw)


openxlsx::write.xlsx(RMLS.only.drivermut.0.1.sig,"RMLS.only.drivermut.0.1.sig.genes.xlsx")

plotEnrich(Enriched_RMLS.only.drivermut.0.1.sig[[1]], showTerms = 30, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis ")


library(openxlsx)

RMLS_only_sig_DP_mutation_genes

write.xlsx(RMLS_only_sig_DP_mutation_genes,"RMLS_only_sig_DP_driver_mutation.xlsx")

RMLS_sig_DP_cytoband_sort <- RMLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
 select(Unique_Name, Cytoband) %>%
 .[grepl( Unique_Name,pattern = "DP",.)]


MLS_sig_DP_cytoband_sort <- MLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband) %>%
  .[grepl( Unique_Name,pattern = "DP",.)]


RMLS_only_sig_DP_cytoband <- anti_join(RMLS_sig_DP_cytoband_sort,MLS_sig_DP_cytoband_sort,"Cytoband")

RMLS_only_sig_DP_cytoband_genes <- RMLS.gistic@data %>% filter(Cytoband %in% RMLS_only_sig_DP_cytoband$Unique_Name) %>% select(Hugo_Symbol, Cytoband) 

RMLS_only_sig_AP_cytoband_genes

RMLS_only_sig_DP_cytoband_genes

RMLS_sig_AP_cytoband_sort <- RMLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband) %>%
  .[grepl( Unique_Name,pattern = "AP",.)]


MLS_sig_AP_cytoband_sort <- MLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband) %>%
  .[grepl( Unique_Name,pattern = "AP",.)]


RMLS_only_sig_AP_cytoband <- anti_join(RMLS_sig_AP_cytoband_sort,MLS_sig_AP_cytoband_sort,"Cytoband")

RMLS_CNV_cytoband <- RMLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband,qvalues) 

MLS_CNV_cytoband <- MLS.gistic@cytoband.summary %>% 
  filter(qvalues<0.25) %>% 
  select(Unique_Name, Cytoband,qvalues)

RMLS_only_sig_all_cytoband <- anti_join(RMLS_CNV_cytoband,MLS_CNV_cytoband,"Cytoband")



RMLS_cn_matrix <- RMLS.gistic@cnMatrix


MLS_cn_matrix <- MLS.gistic@cnMatrix

coOncoplot(m1 = MLS_maf,m2 = RMLS_maf,genes = RMLS_only_sig_DP_mutation_genes$Hugo_Symbol,)



RMLS_only_sig_DP_mutation_genes


library(tibble)

RMLS_cn_matrix_df <- rownames_to_column(data.frame(RMLS_cn_matrix),"Unique_Name")

RMLS_cn_matrix_df <- merge(RMLS_CNV_cytoband,RMLS_cn_matrix_df)


MLS_cn_matrix_df <- rownames_to_column(data.frame(MLS_cn_matrix),"Unique_Name")

MLS_cn_matrix_df <- merge(MLS_CNV_cytoband,MLS_cn_matrix_df)

MLS_cn_matrix_df_for_RMLS <- merge(RMLS_CNV_cytoband,MLS_cn_matrix_df[,-1],
                                   by= "Cytoband",sort = FALSE,all.x =TRUE)


RMLS_cn_only_matrix_df <- merge(RMLS_only_sig_all_cytoband,RMLS_cn_matrix_df)

RMLS_cn_only_matrix_df = RMLS_cn_only_matrix_df[order(RMLS_cn_only_matrix_df$qvalues),]

RMLS_cn_only_matrix_df2 <- RMLS_cn_only_matrix_df[,-c(2:3)]

RMLS_cn_only_matrix_df2 <- column_to_rownames(RMLS_cn_only_matrix_df2,"Unique_Name")


col_cnv <- c("Amp" = "#C11C0F", "Del"="#2166AC")


alter_fun_cnv <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#C11C0F", col = "#C11C0F"))
  },
  # bug red
  Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#2166AC", col = "#2166AC"))
  }
)


heatmap_legend_cnv <- list(title = "", at = c("Amp","Del"),labels = c("Amp","Del"),labels_gp = gpar(fontsize = 5,fontface =2,fontfamily="Arial"),ncol = 2, by_row = TRUE)

library(ggplot2)

library(ggpubr)

library(ComplexHeatmap)

h5 = oncoPrint(RMLS_cn_only_matrix_df2,show_pct = FALSE,right_annotation = gene_anno,
              alter_fun = alter_fun_cnv,col = col_cnv,show_column_names = TRUE,remove_empty_columns = FALSE,remove_empty_rows = FALSE,top_annotation = NULL, heatmap_legend_param = heatmap_legend_cnv, column_names_gp=gpar(fontsize = 5.8,fontface =2,fontfamily="Arial"),row_names_gp=gpar(fontsize = 5.8,fontface =2,fontfamily="Arial"), width = ncol(RMLS_cn_only_matrix_df2)*unit(7, "mm"), height = nrow(RMLS_cn_only_matrix_df2)*unit(3, "mm")) 

#"HLA-A \nTNF"

gene_anno = rowAnnotation(link = anno_mark(at = c(1,2,3,4,5,7,15,16,20,22,25), 
                                   labels = c("TNF","HES5","PTEN","PTPRK","FAN1 \nKLF13","TTN","ZNF460","RAD7","MAP3K4","MYO1C","CDR2"),which= 'row',labels_gp = gpar(fontsize = 7,fontface =2,fontfamily="Arial"), padding = unit(1.5, "lines")))



png(filename="final_result_MUT/RMLS_sig_CNV_plot.png",width=800,height=600,unit="px",bg="transparent",res = 150)

draw(h5,heatmap_legend_side = "bottom")

dev.off()


###############################################################
#oncodrivefml3 conda
#oncodrivefml -i ./MLS_oncodrive_input.txt -e /SSD/jeeh99/tools/oncodrivefml2/resource/cds_hg38_chr.tsv -c /SSD/jeeh99/tools/oncodrivefml2/resource/oncodrivefml_v2.conf --sequencing wes -o ./MLS_re --cores 8


mls.oncodrive.input <- MLS_maf@data[ ,c("Chromosome",
                            "Start_Position",
                            "Reference_Allele",
                            "Tumor_Seq_Allele2",
                            "Tumor_Sample_Barcode"
                            )]


MLS_maf_dt


names(mls.oncodrive.input) <- c("CHROMOSOME", 
                            "POSITION",
                            "REF",
                            "ALT",
                            "SAMPLE")

write.table(mls.oncodrive.input, file="MLS_oncodrive_input.txt", sep = "\t", row.names=FALSE, quote=FALSE)

#mls.oncodrive.input$CHROMOSOME<- gsub("chr","",mls.oncodrive.input$CHROMOSOME)

rmls.oncodrive.input <- RMLS_maf@data[ ,c("Chromosome",
                                        "Start_Position",
                                        "Reference_Allele",
                                        "Tumor_Seq_Allele2",
                                        "Tumor_Sample_Barcode"
)]


names(rmls.oncodrive.input) <- c("CHROMOSOME", 
                                "POSITION",
                                "REF",
                                "ALT",
                                "SAMPLE")

write.table(rmls.oncodrive.input, file="RMLS_oncodrive_input.txt", sep = "\t", row.names=FALSE, quote=FALSE)


MLS.oncodrive.res <- read.table("MLS_oncodrive_input-oncodrivefml.tsv/MLS_oncodrive_input-oncodrivefml.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

MLS.oncodrive.res<- na.omit(MLS.oncodrive.res)

MLS.oncodrive_sig.res =MLS.oncodrive.res[MLS.oncodrive.res$Q_VALUE<0.05,c("SYMBOL","Q_VALUE")]

names(MLS.oncodrive_sig.res) <- c("Hugo_Symbol","oncodrivefml")

MLS.dndscv_sig.res <- signif_genes_MLS

names(MLS.dndscv_sig.res) <- c("Hugo_Symbol","dndscv")

MLS.drivermut.sig <- merge(MLS.dndscv_sig.res,MLS.oncodrive_sig.res,"Hugo_Symbol")
#####################################################################
library(ggplot2)

library(enrichR)

Enriched_go_all_driver_mls <- enrichr(genes = MLS.drivermut.sig$Hugo_Symbol, databases = dbs_go)


png(filename="final_result_MUT/Drivermut/MLS_all_sig_genes_go_bp.png",width=1000,height=800,unit="px",bg="transparent",res=100)


plotEnrich(Enriched_go_all_driver_mls[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for driver mutation in MLS")

dev.off()


Enriched_pw_all_driver_mls <- enrichr(genes = MLS.drivermut.sig$Hugo_Symbol, databases = dbs_pw)


png(filename="final_result_MUT/Drivermut/MLS_all_sig_genes_kegg.png",width=1000,height=800,unit="px",bg="transparent",res=100)

plotEnrich(Enriched_pw_all_driver_mls[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for driver mutation in MLS")

dev.off()


Enriched_ppi_all_driver_mls <- enrichr(genes =  MLS.drivermut.sig$Hugo_Symbol, databases = dbs_ppi)


png(filename="final_result_MUT/Drivermut/MLS_all_sig_genes_ppi.png",width=1000,height=800,unit="px",bg="transparent",res=100)

plotEnrich(Enriched_ppi_all_driver_mls[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "PPI enrichment analysis for driver mutation in MLS")

dev.off()

###############################################################

RMLS.oncodrive.res <- read.table("RMLS_oncodrive_input-oncodrivefml.tsv/RMLS_oncodrive_input-oncodrivefml.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

RMLS.oncodrive.res<- na.omit(RMLS.oncodrive.res)

RMLS.oncodrive_sig.res =RMLS.oncodrive.res[RMLS.oncodrive.res$Q_VALUE<0.05,c("SYMBOL","Q_VALUE")]

names(RMLS.oncodrive_sig.res) <- c("Hugo_Symbol","oncodrivefml")

RMLS.dndscv_sig.res <- RMLS_signif_genes 

names(RMLS.dndscv_sig.res) <- c("Hugo_Symbol","dndscv")

RMLS.drivermut.sig <- merge(RMLS.dndscv_sig.res,RMLS.oncodrive_sig.res,"Hugo_Symbol")
#################################################################

Enriched_go_all_driver_rmls <- enrichr(genes = RMLS.drivermut.sig$Hugo_Symbol, databases = dbs_go)


png(filename="final_result_MUT/Drivermut/RMLS_all_sig_genes_go_bp.png",width=1000,height=800,unit="px",bg="transparent",res=100)


plotEnrich(Enriched_go_all_driver_rmls[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for driver mutation in RMLS")

dev.off()


Enriched_pw_all_driver_rmls <- enrichr(genes = RMLS.drivermut.sig$Hugo_Symbol, databases = dbs_pw)

library(ggthemes)
#,width=800,height=600,unit="px",bg="transparent",res = 150
#width=1000,height=800,unit="px",bg="transparent",res=100

library(enrichplot)
library(enrichR)

png(filename="final_result_MUT/Drivermut/RMLS_all_sig_genes_kegg2.png",width=800,height=600,unit="px",bg="transparent",res = 100)

plotEnrich(Enriched_pw_all_driver_rmls[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for driver mutation in RMLS") + theme_few()+
  theme(axis.text.x = element_text(face = "bold",size =8,family = "arial",colour = "black"),                      # adjusting the position
        axis.title.y = element_text(face = "bold",size = 8,family = "arial"),                   # face the x axit title/label
        axis.text.y = element_text(face = "bold",size = 8,family = "arial",colour = "black"),                   # face the y axis title/label
        plot.title = element_text(hjust = 0.1,face = "bold",size =10 )) + scale_fill_gradient(low ="#C11C0F", high = "#2166AC") + ylab("Gene Count")+
  xlab("")


dev.off()

KEGG_driver_mutation_RMLS <- Enriched_pw_all_driver_rmls[[1]]

Enriched_ppi_all_driver_rmls <- enrichr(genes =  RMLS.drivermut.sig$Hugo_Symbol, databases = dbs_ppi)


png(filename="final_result_MUT/Drivermut/RMLS_all_sig_genes_ppi.png",width=1000,height=800,unit="px",bg="transparent",res=100)

dd = plotEnrich(Enriched_ppi_all_driver_rmls[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")

dd$layers

dev.off()



KEGG_driver_mutation_RMLS[, .(count=.N), by = Genes]

KEGG_driver_mutation_RMLS %>%
  dplyr:: count(Genes, name = 'Gene_count')

library(dplyr)

top20_MLS_vaf_re <- merge(top20_mut_df,top20_MLS_vaf,sort = FALSE,all.x = TRUE)


MLS_all_res_genes_dndscv <-  MLS_sel_cv[,c("gene_name","qglobal_cv")]

names(MLS_all_res_genes_dndscv) <- c("Hugo_Symbol","dndscv")

MLS_all_res_genes_dndscv$dndscv <- -log10(MLS_all_res_genes_dndscv$dndscv)


top20_MLS_vaf_re <- merge(top20_MLS_vaf,MLS_all_res_genes_dndscv,sort = FALSE,all.x = TRUE)


MLS_all_res_genes_oncodrivefml <-  MLS.oncodrive.res[,c("SYMBOL","Q_VALUE")]


names(MLS_all_res_genes_oncodrivefml) <- c("Hugo_Symbol","oncodrivefml")


MLS_all_res_genes_oncodrivefml$oncodrivefml <- -log10(MLS_all_res_genes_oncodrivefml$oncodrivefml)


top20_MLS_vaf_re <- merge(top20_MLS_vaf_re ,MLS_all_res_genes_oncodrivefml,sort = FALSE,all.x = TRUE)


top20_MLS_vaf_re[is.na(top20_MLS_vaf_re)] = 0

#######################################################################
top20_RMLS_vaf_re <- merge(top20_mut_df,top20_RMLS_vaf,sort = FALSE,all.x = TRUE)


RMLS_all_res_genes_dndscv <-  RMLS_sel_cv[,c("gene_name","qglobal_cv")]

names(RMLS_all_res_genes_dndscv) <- c("Hugo_Symbol","dndscv")

RMLS_all_res_genes_dndscv$dndscv <- -log10(RMLS_all_res_genes_dndscv$dndscv)


top20_RMLS_vaf_re <- merge(top20_RMLS_vaf,RMLS_all_res_genes_dndscv,sort = FALSE,all.x = TRUE)


RMLS_all_res_genes_oncodrivefml <-  RMLS.oncodrive.res[,c("SYMBOL","Q_VALUE")]


names(RMLS_all_res_genes_oncodrivefml) <- c("Hugo_Symbol","oncodrivefml")


RMLS_all_res_genes_oncodrivefml$oncodrivefml <- -log10(RMLS_all_res_genes_oncodrivefml$oncodrivefml)


top20_RMLS_vaf_re <- merge(top20_RMLS_vaf_re ,RMLS_all_res_genes_oncodrivefml,sort = FALSE,all.x = TRUE)


top20_RMLS_vaf_re[is.na(top20_RMLS_vaf_re)] = 0


top20_RMLS_vaf_re <- merge(top20_mut_df,top20_RMLS_vaf_re,sort = FALSE,all.x = TRUE)
###############################################################################
library(ComplexHeatmap)

col <- c("Nonsense_Mutation" = "red","Missense_Mutation" = "#008000","In_Frame_Del"="yellow","Frame_Shift_Del"="blue","Translation_Start_Site" ="orange","Multi_Hit"="black")


col2 <- c("Nonsense_Mutation" = "red","Missense_Mutation" = "#008000","In_Frame_Del"="yellow","Frame_Shift_Del"="blue","Multi_Hit"="black")


alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  }
)



alter_fun2 <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  # small green
  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),  h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Multi_Hit"], col = NA))
  }
)


heatmap_legend_param <- list(title = "", at = c("Nonsense_Mutation" , "Missense_Mutation","In_Frame_Del","Frame_Shift_Del","Translation_Start_Site","Multi_Hit"),labels = c("Nonsense_Mutation","Missense_Mutation", 
                                                                                                                                                                            "In_Frame_Del","Frame_Shift_Del","Translation_Start_Site","Multi_Hit"),labels_gp = gpar(fontsize = 5),ncol = 3, by_row = TRUE)



heatmap_legend_param2 <- list(title = "", at = c("Nonsense_Mutation" , "Missense_Mutation","In_Frame_Del",
                                                 "Frame_Shift_Del","Multi_Hit"), 
                              labels = c("Nonsense_Mutation","Missense_Mutation", 
                                         "In_Frame_Del","Frame_Shift_Del","Multi_Hit"))


h1 = rowAnnotation("oncodrivefml" = anno_barplot(top20_MLS_vaf_re$oncodrivefml,axis_param=list(direction = "reverse",gp=gpar(fontsize = 4),at = c(0, 1, 2))),width = unit(1.0, "cm"), annotation_name_gp = gpar(fontsize=5)) + rowAnnotation("dndscv" = anno_barplot(top20_MLS_vaf_re$dndscv,axis_param=list(direction = "reverse",gp=gpar(fontsize = 4),at = c(0, 1, 2))),width = unit(1.0, "cm"), annotation_name_gp = gpar(fontsize=5)) + oncoPrint(MLS_mut_mat ,show_pct = FALSE,right_annotation =rowAnnotation(Type= c("Both",rep("MLS",5),"Both",rep("MLS",5),rep("RMLS",10)), col = list(Type = c("MLS" = "blue","RMLS" = "red" ,"Both" = "green")),annotation_legend_param = list(Type = list(title = "",labels_gp = gpar(fontsize = 5),ncol = 3, by_row = TRUE)), width = unit(1, "mm"), annotation_name_gp = gpar(fontsize=5)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                       left_annotation = rowAnnotation("VAF" = anno_barplot(top20_MLS_vaf$VAF,axis_param = list(direction = "reverse",gp=gpar(fontsize = 4),at = c(0, 0.2, 0.4))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       width = unit(1.0, "cm"),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       annotation_name_gp = gpar(fontsize=5)),
                                                                                                                                                                                                                                                                                                                                                                                                                                                       alter_fun = alter_fun, col = col,
                                                                                                                                                                                                                                                                                                                                                                                                                                                       column_labels = colnames(MLS_mut_mat),
                                                                                                                                                                                                                                                                                                                                                                                                                                                       show_column_names = TRUE,
                                                                                                                                                                                                                                                                                                                                                                                                                                                       remove_empty_columns = FALSE,
                                                                                                                                                                                                                                                                                                                                                                                                                                                       remove_empty_rows = FALSE,top_annotation = NULL,
                                                                                                                                                                                                                                                                                                                                                                                                                                                       heatmap_legend_param = heatmap_legend_param,
                                                                                                                                                                                                                                                                                                                                                                                                                                                       column_names_gp=gpar(fontsize = 5),
                                                                                                                                                                                                                                                                                                                                                                                                                                                       row_names_gp=gpar(fontsize = 5))  



h2 = oncoPrint(RMLS_mut_mat,show_row_names = FALSE, show_pct = FALSE, right_annotation = rowAnnotation("VAF" = anno_barplot(top20_RMLS_vaf$VAF,axis_param=list(gp=gpar(fontsize = 4),at = c(0, 0.2, 0.4)),width = unit(1.0, "cm")),                                                             annotation_name_gp = gpar(fontsize=5)),
               alter_fun = alter_fun2, col = col2, column_labels = colnames(RMLS_mut_mat) ,show_column_names = TRUE,
               left_annotation = NULL,
               remove_empty_columns = FALSE, top_annotation = NULL,
               remove_empty_rows = FALSE ,row_labels = NULL,show_heatmap_legend = FALSE,column_names_gp=gpar(fontsize = 5)) + rowAnnotation("dndscv" = anno_barplot(top20_RMLS_vaf_re$dndscv,axis_param=list(gp=gpar(fontsize = 4),at = c(0, 1, 2))),width = unit(1.0, "cm"), annotation_name_gp = gpar(fontsize=5)) +
  rowAnnotation("oncodrivefml" = anno_barplot(top20_RMLS_vaf_re$oncodrivefml,axis_param=list(gp=gpar(fontsize = 4),at = c(0, 1, 2))),width = unit(1.0, "cm"), annotation_name_gp = gpar(fontsize=5))

png(filename="final_result_MUT/Drivermut/top20_coonco_plot.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

draw(h1+ h2,heatmap_legend_side = "bottom")

dev.off()


#####################################################################

library(dplyr)
library(openxlsx)

RMLS.only.drivermut.sig <- anti_join(RMLS.drivermut.sig,MLS.drivermut.sig,"Hugo_Symbol")

write.xlsx(RMLS.only.drivermut.0.1.sig,"RMLS.only.drivermut.0.1.xlsx")

RMLS.only.drivermut.sig

Enriched_only_rmls_dr_go <- enrichr(genes = RMLS.only.drivermut.sig$Hugo_Symbol, databases = dbs_go)


png(filename="final_result_MUT/Drivermut/RMLS_only_sig_genes_go_bp.png",width=1000,height=800,unit="px",bg="transparent",res=100)


plotEnrich(Enriched_only_rmls_dr_go[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "GO BP enrichment analysis for driver mutation in only RMLS")

dev.off()

Enriched_only_rmls_dr_pw <- enrichr(genes = RMLS.only.drivermut.sig$Hugo_Symbol, databases = dbs_pw)


png(filename="final_result_MUT/Drivermut/RMLS_only_sig_genes_kegg.png",width=1000,height=800,unit="px",bg="transparent",res=100)


plotEnrich(Enriched_only_rmls_dr_pw [[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "KEGG enrichment analysis for driver mutation in only RMLS")

dev.off()


png(filename="final_result_MUT/Drivermut/RMLS_only_sig_genes_wiki.png",width=1000,height=800,unit="px",bg="transparent",res=100)


plotEnrich(Enriched_only_rmls_dr_pw [[2]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "WiKipathways enrichment analysis for driver mutation in only RMLS")

dev.off()

Enriched_only_rmls_dr_ppi <- enrichr(genes = RMLS.only.drivermut.sig$Hugo_Symbol, databases = dbs_ppi)

png(filename="final_result_MUT/Drivermut/RMLS_only_sig_genes_ppi.png",width=1000,height=800,unit="px",bg="transparent",res=100)


plotEnrich(Enriched_only_rmls_dr_ppi[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value",title = "PPI enrichment analysis for driver mutation in only RMLS")

dev.off()

library(openxlsx)

write.xlsx(RMLS.only.drivermut.sig,"RMLS.only.drivermutation_genes.xlsx")

write.xlsx(RMLS.drivermut.sig,"RMLS.sig.drivermutation_genes.xlsx") 
##################################################################
MLS.dndscv_sig.0.1.res <- MLS_sel_cv[MLS_sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

names(MLS.dndscv_sig.0.1.res) <- c("Hugo_Symbol","dndscv")

MLS.oncodrive_sig.0.1.res =MLS.oncodrive.res[MLS.oncodrive.res$Q_VALUE<0.1,c("SYMBOL","Q_VALUE")]

names(MLS.oncodrive_sig.0.1.res) <- c("Hugo_Symbol","oncodrivefml")

MLS.drivermut.0.1.sig <- merge(MLS.dndscv_sig.0.1.res,
                               MLS.oncodrive_sig.0.1.res
                               ,"Hugo_Symbol")

RMLS.dndscv_sig.0.1.res <- RMLS_sel_cv[RMLS_sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

names(RMLS.dndscv_sig.0.1.res) <- c("Hugo_Symbol","dndscv")

RMLS.oncodrive_sig.0.1.res =RMLS.oncodrive.res[RMLS.oncodrive.res$Q_VALUE<0.1,c("SYMBOL","Q_VALUE")]

names(RMLS.oncodrive_sig.0.1.res) <- c("Hugo_Symbol","oncodrivefml")

RMLS.drivermut.0.1.sig <- merge(RMLS.dndscv_sig.0.1.res,
                               RMLS.oncodrive_sig.0.1.res
                               ,"Hugo_Symbol")


library(dplyr)

RMLS.only.drivermut.0.1.sig <- anti_join(RMLS.drivermut.0.1.sig,MLS.drivermut.0.1.sig,"Hugo_Symbol")



MLS.only.drivermut.0.1.sig <- anti_join(MLS.drivermut.0.1.sig,RMLS.drivermut.0.1.sig,"Hugo_Symbol")


library(openxlsx)

write.xlsx(RMLS.drivermut.0.1.sig,"RMLS.sig.0.1.drivermutation_genes.xlsx")



write.xlsx(MLS.drivermut.0.1.sig,"MLS.sig.0.1.drivermutation_genes.xlsx")
library(maftools)

lollipopPlot(
  maf = RMLS_all_maf,
  gene = 'PTEN',
  AACol = 'aaChange',
  showMutationRate = TRUE
)


lollipopPlot(
  maf = RMLS_maf,
  gene = 'MICAL3',
  AACol = 'aaChange',
  showMutationRate = TRUE
)


lollipopPlot(
  maf = MLS_maf,
  gene = 'PIK3CA',
  AACol = 'aaChange',
  showMutationRate = TRUE
)


MLS_all_maf_df <- MLS_all_maf@data

RMLS_all_maf_df <- RMLS_all_maf@data


library(fgsea)

KEGG_pathway<- fgsea::gmtPathways("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")


KEGG_pathway$KEGG_HEDGEHOG_SIGNALING_PATHWAY

KEGG_pathway$KEGG_CELL_CYCLE

KEGG_pathway$KEGG_MTOR_SIGNALING_PATHWAY

KEGG_pathway$KEGG_BASE_EXCISION_REPAIR

KEGG_pathway$KEGG_MISMATCH_REPAIR

KEGG_pathway$KEGG_NUCLEOTIDE_EXCISION_REPAIR

KEGG_pathway$KEGG_FATTY_ACID_METABOLISM

result_pathway = data.frame(Gene = c(KEGG_pathway$KEGG_HEDGEHOG_SIGNALING_PATHWAY,KEGG_pathway$KEGG_CELL_CYCLE,KEGG_pathway$KEGG_MTOR_SIGNALING_PATHWAY,KEGG_pathway$KEGG_BASE_EXCISION_REPAIR,KEGG_pathway$KEGG_MISMATCH_REPAIR,KEGG_pathway$KEGG_NUCLEOTIDE_EXCISION_REPAIR,KEGG_pathway$KEGG_FATTY_ACID_METABOLISM
                                     ),
                            Pathway = rep(c(
  "Hedgehog Signaling", "Cell Cycle", "MTOR Signalling", "BER", "MMR","NER","Fatty Acid Metabolism"
), c(56, 125, 52, 35, 23,44,42)),
stringsAsFactors = FALSE)

oncoplot(maf = RMLS_maf, pathways = "sigpw" , gene_mar = 8, fontSize = 0.6, topPathways = 5, collapsePathway = TRUE,showTumorSampleBarcodes = T)


oncoplot(maf = MLS_maf, pathways = "sigpw" , gene_mar = 8, fontSize = 0.6, topPathways = 5, collapsePathway = TRUE,showTumorSampleBarcodes = T)

###########################################################################


RMLS_only_sig_DP_mutation_genes


oncoplot(maf = RMLS_all_maf,genes = RMLS_only_sig_DP_mutation_genes
$Hugo_Symbol,keepGeneOrder = TRUE,writeMatrix = T)


RMLS_mut_sig_cnv <- read.table("onco_matrix.txt",sep = "\t",stringsAsFactors = FALSE,header = TRUE)



col7 <- c("Nonsense_Mutation" = "red","Missense_Mutation" = "#008000","In_Frame_Del"="yellow","Frame_Shift_Del"="#CC79A7","Multi_Hit"="black","Del"="#2166AC")



alter_fun7 <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.9, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.9, "mm"), 
              gp = gpar(fill = col7["Nonsense_Mutation"], col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.9, "mm"), 
              gp = gpar(fill = col7["Missense_Mutation"], col = NA))
  },
  # small green
  
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.9, "mm"), 
              gp = gpar(fill = col7["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(0.9, "mm"), 
              gp = gpar(fill = col7["Frame_Shift_Del"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"),  h-unit(0.9, "mm"), 
              gp = gpar(fill = col7["Multi_Hit"], col = NA))
  },
  Del = function(x, y, w, h) {
    grid.polygon(
      unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
      unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
      gp = gpar(fill = col7["Del"], col = "white"))
  }
)


heatmap_legend_param7 <- list(title = "", at = c("Nonsense_Mutation" , "Missense_Mutation","In_Frame_Del","Frame_Shift_Del","Multi_Hit","Del"),labels = c("Nonsense_Mutation","Missense_Mutation",  "In_Frame_Del","Frame_Shift_Del","Multi_Hit","Del"),labels_gp = gpar(fontsize = 5.8,fontface = 2,fontfamily="Arial"),ncol = 3, by_row = TRUE)

h7 <- oncoPrint(RMLS_mut_sig_cnv,show_pct = FALSE,
                alter_fun = alter_fun7,col = col7,show_column_names = TRUE,remove_empty_columns = FALSE,remove_empty_rows = FALSE,top_annotation = NULL, heatmap_legend_param =heatmap_legend_param7, column_names_gp=gpar(fontsize = 5.8,fontface =2,fontfamily="Arial"),row_names_gp=gpar(fontsize = 5.8,fontface =2,fontfamily="Arial")) 


png(filename="final_result_MUT/RMLS_onco_cnv_plot.png",width=800,height=600,unit="px",bg="transparent",res = 150)


draw(h7,heatmap_legend_side = "bottom")

dev.off()
