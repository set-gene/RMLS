library(BSgenome)
library(MutationalPatterns)

library(ggplot2)

library(RColorBrewer)

names <- c ("MLS_1","MLS_2","MLS_3", "MLS_4", "MLS_5", "MLS_6","MLS_7","MLS_8","MLS_9",
            "RMLS_1","RMLS_2","RMLS_3", "RMLS_4", "RMLS_5", "RMLS_6","RMLS_7","RMLS_8","RMLS_9")

files <-list.files(path = "./",pattern= ".vcf", full.names = TRUE)

files

ref_genome <-"BSgenome.Hsapiens.UCSC.hg38"

library(ref_genome, character.only = TRUE)

rmmls <- read_vcfs_as_granges(files,names,ref_genome)

muts = mutations_from_vcf(rmmls[[1]])

type_occurrences = mut_type_occurrences(rmmls, ref_genome)

plot_spectrum(type_occurrences ,CT=TRUE,legend = TRUE)

write.table(type_occurrences,"ALL_SNP.txt", sep = "\t")

tissue <- c ("MLS_1","MLS_2","MLS_3", "MLS_4", "MLS_5", "MLS_6","MLS_7","MLS_8","MLS_9",
                 "RMLS_1","RMLS_2","RMLS_3", "RMLS_4", "RMLS_5", "RMLS_6","RMLS_7","RMLS_8","RMLS_9")

tissue_group =  c(rep("MLS", 9), rep("RMLS", 9))


signature_group =  c(rep("A", 5), rep("B", 4),rep("A", 5), rep("B", 4))


signature_tissue_group =  c(rep("MLS_A", 5), rep("MLS_B", 4),rep("RMLS_A", 5), rep("RMLS_B", 4))


mut_mat <- mut_matrix(vcf_list = rmmls, ref_genome = ref_genome)

nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 10)


nmf_res2 <- extract_signatures(mut_mat, rank = 4, nrun = 10)

colnames(nmf_res2$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")


rownames(nmf_res2$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")


colnames(nmf_res$signatures) <- c("Signature A", "Signature B")

rownames(nmf_res$contribution) <- c("Signature A", "Signature B")


png(filename="final_result_MUT/all_signature_AB.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

plot_96_profile(nmf_res$signatures, condensed = TRUE)

dev.off()

plot_96_profile(nmf_res2$signatures, condensed = TRUE)


library(gridExtra)

pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")

pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "absolute")

grid.arrange(pc1, pc2)

png(filename="final_result_MUT/all_signature.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative", coord_flip = TRUE,palette = c("skyblue","red"))

dev.off()

pch1 <-plot_contribution_heatmap(nmf_res$contribution,sig_order = c("Signature B","Signature A"))


pch2 <-plot_contribution_heatmap(nmf_res2$contribution,sig_order = c("Signature A","Signature B","Signature C","Signature D"))


pch1 <-plot_contribution_heatmap(nmf_res$contribution,sig_order = c("Signature B","Signature A"))

grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6))

library(MutationalPatterns)

png(filename="final_result_MUT/all_of_signature.png",width=2000,height=1500,unit="px",bg="transparent",res = 100)

plot_spectrum(type_occurrences,by = tissue , CT = TRUE, legend = TRUE)

dev.off()

png(filename="final_result_MUT/each_signature.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plot_spectrum(type_occurrences, by = signature_tissue_group , CT = TRUE, legend = TRUE)

dev.off()


png(filename="final_result_MUT/group_signature.png",width=800,height=600,unit="px",bg="transparent",res = 100)


plot_spectrum(type_occurrences, by = tissue_group , CT = FALSE, legend = TRUE,colors =  c( "#C11C0F","#1b9e77","#2166AC","#CC79A7","#D55E00","#e6ab02","grey50" )) + theme(axis.text.x = element_text(face = "bold",size = 10,family = "arial",colour = "black"), axis.title.y = element_text(face = "bold",size = 11,family = "arial"),axis.text.y = element_text(face = "bold",size = 12,family = "arial",colour = "black"), axis.line = element_line(size=0.5, colour = "black"),legend.text = element_text(face = "bold",size = 11,family = "arial"),legend.title = element_text(face = "bold",size = 12,family = "arial"),text= element_text(face = "bold",size = 12,family = "arial"))


dev.off()


grid.arrange(p4, ncol = 2)


sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")

cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)

cancer_signatures$Somatic.Mutation.Type

new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)

cancer_signatures$Somatic.Mutation.Type

cancer_signatures = cancer_signatures[as.vector(new_order),]

row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type

cancer_signatures = as.matrix(cancer_signatures[,4:33])

cancer_signatures


hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")

cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]

cosmic_order 

plot(hclust_cosmic)

cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)

cos_sim_samples_signatures_t = t(cos_sim_samples_signatures)

cos_sim_samples_signatures_t

plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE)

plot_cosine_heatmap(cos_sim_samples_signatures_t,col_order = cosmic_order,cluster_rows = TRUE)

cos_sim_samples_signatures


fit_res <-fit_to_signatures(mut_mat, cancer_signatures)

select <- which(rowSums(fit_res$contribution) > 10)


png(filename="final_result_MUT/all_cosmic_signature.png",width=1500,height=1000,unit="px",bg="transparent",res = 100)


plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = TRUE,mode = "relative")

dev.off()


select_sig = subset(select,select == 1 | select == 4 | select == 6 | select == 20 | select == 24)

select_sig

select_sig = data.frame(select)

select_sig = t(select_sig)

select_sig[,29]


mypalette <- colorRampPalette(brewer.pal(8,"Set3"))

colorRampPalette(brewer.pal())

mypalette2 <- colorRampPalette(brewer.pal(8,"Set2"))

plot_contribution(fit_res$contribution[select_sig ,],cancer_signatures[,select_sig ],coord_flip = TRUE,mode = "relative")

fit_res$contribution

cancer_signatures


df = fit_res$contribution[, apply(fit_res$contribution, 2, function(x){sum(x>0)})> 5]


############################################################2 signature group contribution to COSMIC

nmf_res$signatures


###########################################################group signature
pooled_mut_mat <- pool_mut_mat(mut_mat, grouping = tissue_group)

pooled_mut_mat

group_nmf_res <- extract_signatures(pooled_mut_mat, rank = 2, nrun = 10)

group_nmf_res

colnames(group_nmf_res$signatures) <- c("Signature A", "Signature B")

rownames(group_nmf_res$contribution) <- c("Signature A", "Signature B")


png(filename="final_result_MUT/group_signature_AB.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

plot_96_profile(group_nmf_res$signatures, condensed = TRUE,colors =  c( "#C11C0F","#1b9e77","#2166AC","#CC79A7","#D55E00","#e6ab02" ))+ theme(axis.text.x = element_text(face = "bold",size = 5,family = "arial",colour = "black"),axis.title.x = element_text(face = "bold",size = 11,family = "arial") ,axis.title.y = element_text(face = "bold",size = 11,family = "arial"),axis.text.y = element_text(face = "bold",size = 12,family = "arial",colour = "black") ,axis.line = element_line(size=0.5, colour = "black"),legend.text = element_text(face = "bold",size = 11,family = "arial"),legend.title = element_text(face = "bold",size = 12,family = "arial"),text= element_text(face = "bold",size = 12,family = "arial"))

dev.off()

library(MutationalPatterns)

library(ggplot2)

png(filename="final_result_MUT/group_signature_contribution.png",width=800,height=600,unit="px",bg="transparent",res = 100)


plot_contribution(group_nmf_res$contribution, group_nmf_res$signatures,mode = "relative", coord_flip = TRUE,palette = c( "#C11C0F","#2166AC")) + theme(axis.text.x = element_text(face = "bold",size = 10,family = "arial",colour = "black"), axis.title.x = element_text(face = "bold",size = 11,family = "arial"),axis.text.y = element_text(face = "bold",size = 12,family = "arial",colour = "black"), axis.line = element_line(size=0.5, colour = "black"),legend.text = element_text(face = "bold",size = 11,family = "arial"),legend.title = element_text(face = "bold",size = 12,family = "arial"))

dev.off()
#########

cancer_signatures_FG = read.table(sp_url, sep = "\t", header = TRUE)

cancer_signatures_FG

group_order = match(row.names(pooled_mut_mat), cancer_signatures_FG$Somatic.Mutation.Type)

group_order

cancer_signatures_FG = cancer_signatures_FG[as.vector(group_order),]

row.names(cancer_signatures_FG) = cancer_signatures_FG$Somatic.Mutation.Type

cancer_signatures_FG

cancer_signatures_FG = as.matrix(cancer_signatures_FG[,4:33])


hclust_cosmic_FG = cluster_signatures(cancer_signatures_FG, method = "average")

plot(hclust_cosmic_FG)


cosmic_order_FG = colnames(cancer_signatures_FG)[hclust_cosmic_FG$order]

cosmic_order_FG


cos_sim_groups_signatures = cos_sim_matrix(pooled_mut_mat, cancer_signatures_FG)


plot_cosine_heatmap(cos_sim_groups_signatures,cluster_rows = TRUE,cluster_cols = FALSE)

cos_sim_groups_signatures

cos_sim_groups_signatures_t = t(cos_sim_groups_signatures)

cos_sim_groups_signatures_t


fit_res_FG <-fit_to_signatures(pooled_mut_mat, cancer_signatures_FG)


select_FG <- which(rowSums(fit_res_FG$contribution) > 10)


png(filename="final_result_MUT/group_cosmic_signature.png",width=1000,height=800,unit="px",bg="transparent",res = 100)


plot_contribution(fit_res_FG$contribution[select_FG,],cancer_signatures_FG[,select_FG],coord_flip = TRUE,mode = "relative") + scale_fill_brewer(palette="Paired")  + theme(axis.text.x = element_text(face = "bold",size = 10,family = "arial",colour = "black"), axis.title.x = element_text(face = "bold",size = 11,family = "arial"),axis.text.y = element_text(face = "bold",size = 12,family = "arial",colour = "black"), axis.line = element_line(size=0.5, colour = "black"),legend.text = element_text(face = "bold",size = 11,family = "arial"),legend.title = element_text(face = "bold",size = 12,family = "arial"))


dev.off()

#############################################################################
###########################################################signature_tissue_group
pooled_mut_mat2 <- pool_mut_mat(mut_mat, grouping = signature_tissue_group)

pooled_mut_mat2

group_nmf_res2 <- extract_signatures(pooled_mut_mat2, rank = 4, nrun = 10)

group_nmf_res2

colnames(group_nmf_res2$signatures) <- c("Signature A", "Signature B","Signature C","Signature D")

rownames(group_nmf_res2$contribution) <- c("Signature A", "Signature B","Signature C", "Signature D")


png(filename="final_result_MUT/group_signature_ABCD.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

plot_96_profile(group_nmf_res2$signatures, condensed = TRUE)

dev.off()

png(filename="final_result_MUT/subgroup_signature_contribution.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

plot_contribution(group_nmf_res2$contribution, group_nmf_res2$signatures,mode = "relative", coord_flip = TRUE,palette = c("skyblue","red","yellow","green"))

dev.off()
#########

cancer_signatures_FG2 = read.table(sp_url, sep = "\t", header = TRUE)

cancer_signatures_FG2

group_order2 = match(row.names(pooled_mut_mat2), cancer_signatures_FG2$Somatic.Mutation.Type)


cancer_signatures_FG2 = cancer_signatures_FG2[as.vector(group_order2),]

row.names(cancer_signatures_FG2) = cancer_signatures_FG2$Somatic.Mutation.Type

cancer_signatures_FG2

cancer_signatures_FG2 = as.matrix(cancer_signatures_FG2[,4:33])


hclust_cosmic_FG2 = cluster_signatures(cancer_signatures_FG2, method = "average")

plot(hclust_cosmic_FG)


cosmic_order_FG2 = colnames(cancer_signatures_FG2)[hclust_cosmic_FG2$order]

cosmic_order_FG


cos_sim_groups_signatures2 = cos_sim_matrix(pooled_mut_mat2, cancer_signatures_FG2)


plot_cosine_heatmap(cos_sim_groups_signatures2,cluster_rows = TRUE,cluster_cols = FALSE)

cos_sim_groups_signatures

cos_sim_groups_signatures_t = t(cos_sim_groups_signatures)

cos_sim_groups_signatures_t


fit_res_FG2 <-fit_to_signatures(pooled_mut_mat2, cancer_signatures_FG2)


select_FG2 <- which(rowSums(fit_res_FG2$contribution) > 10)


png(filename="final_result_MUT/sub_group_cosmic_signature.png",width=1500,height=1000,unit="px",bg="transparent",res = 100)


plot_contribution(fit_res_FG2$contribution[select_FG2,],cancer_signatures_FG2[,select_FG2],coord_flip = TRUE,mode = "relative")


dev.off()




dat_A[,select]
