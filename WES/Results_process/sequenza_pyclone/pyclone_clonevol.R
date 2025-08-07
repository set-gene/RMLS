library(clonevol)

library(tidyverse)

library(dplyr)

library(openxlsx)

library(data.table)

RMLS_sig_driver_genes <- read.xlsx("RMLS.sig.0.1.drivermutation_genes.xlsx")

pyclone.loci<-read.table("Patient_1_loci.tsv",sep="\t", header=T, stringsAsFactors = F)

cluster = filter(as.data.table(table(pyclone.loci$cluster_id)), N > 2)


pyclone.file <- filter(pyclone.loci,pyclone.loci$cluster_id %in% cluster$V1)

tmp=pyclone.file[,c("mutation_id","cluster_id")]

tmp=tmp[!duplicated(tmp),]

rownames(tmp)=tmp$mutation_id

cluster_clonevol=matrix(NA,nrow=length(unique(pyclone.file$mutation_id)),ncol=length(unique(pyclone.file$sample_id)))

rownames(cluster_clonevol)=unique(pyclone.file$mutation_id)

colnames(cluster_clonevol)=unique(pyclone.file$sample_id)

for (i in 1:dim(pyclone.file)[1])
{cluster_clonevol[pyclone.file$mutation_id[i],pyclone.file$sample_id[i]]=pyclone.file$cellular_prevalence[i]*100}

cluster_clonevol=data.frame(cluster=tmp[rownames(cluster_clonevol),"cluster_id"],
                            cluster_clonevol,stringsAsFactors=F)


cluster_clonevol$cluster <- cluster_clonevol$cluster + 1

cluster_clonevol <- cluster_clonevol[,c("cluster","MLS_1","RMLS_1")]

sample.groups = c("MLS","RMLS")
  
names(sample.groups) = colnames( cluster_clonevol)[-1]

sample.groups

library(tidyr)

##################################################

y = infer.clonal.models(variants = cluster_clonevol, cluster.col.name = 'cluster' ,
                        ccf.col.names = colnames( cluster_clonevol)[-1],cancer.initiation.model='polyclonal',
                        subclonal.test = 'bootstrap',subclonal.test.model = 'non-parametric',sample.groups = sample.groups,
                        num.boots = 1000,cluster.center = 'mean',
                        clone.colors = NULL,min.cluster.vaf = 0.01,sum.p = 0.05,alpha = 0.05,random.seed=1)


###################################################################################################
library(tibble)

cluster_clonevol_df <- rownames_to_column(cluster_clonevol,"mutation_id")


library("splitstackshape")

cluster_clonevol_df <- cSplit(indt = cluster_clonevol_df, splitCols = "mutation_id", sep = ":",drop = FALSE) 

colnames(cluster_clonevol_df)[5:7] <- c("CHROM","POS","Hugo_Symbol")

cluster_clonevol_df$mutation_id <- cluster_clonevol_df$Hugo_Symbol

cluster_clonevol_df<- cluster_clonevol_df[,-c(5:7)]



#######################################################
###########################################

cluster_clonevol_df$driver <- ifelse(cluster_clonevol_df$mutation_id %in% RMLS_sig_driver_genes$Hugo_Symbol,TRUE,FALSE)

clone.tree<- transfer.events.to.consensus.trees(y, cluster_clonevol_df[cluster_clonevol_df$driver,],cluster.col.name = 'cluster', event.col.name = 'mutation_id')

clone.tree.branch <- convert.consensus.tree.clone.to.branch(clone.tree, branch.scale = "sqrt")

plot.clonal.models(clone.tree.branch,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'mutation_id',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 0.5,
                   bell.curve.step = 1,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 20,
                   tree.node.text.size = 1,
                   merged.tree.node.size.scale = 1,
                   merged.tree.node.text.size.scale =1,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   ####################tree gene text size
                   mtcab.branch.text.size = 0.4,
                   mtcab.branch.width = 0.4,
                   mtcab.node.size = 2,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 0.6,
                   mtcab.branch.angle = 25,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = './Patient_1',
                   overwrite.output = TRUE,
                   width = 25,
                   height = 10,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(6,8,6,8,8))

##################################################################################################
  library(tibble)
  
  cluster_clonevol_df2<- rownames_to_column(cluster_clonevol,"mutation_id")
  
  
  library("splitstackshape")
  
  cluster_clonevol_df2 <- cSplit(indt = cluster_clonevol_df2, splitCols = "mutation_id", sep = ":",drop = FALSE) 
  
  colnames(cluster_clonevol_df2)[5:7] <- c("CHROM","POS","Hugo_Symbol")
  
  cluster_clonevol_df2$mutation_id <- cluster_clonevol_df2$Hugo_Symbol
  
  
  cluster_clonevol_df2$driver <- ifelse(cluster_clonevol_df2$mutation_id %in% RMLS_sig_driver_genes$Hugo_Symbol,TRUE,FALSE)

cluster_clonevol_df2 <- cluster_clonevol_df2[,c("CHROM","POS","Hugo_Symbol","MLS_1","RMLS_1","driver","cluster")]



write.table(driver_pyclone,"pyclone_driver_genes_p1.txt",sep = "\t",quote = FALSE,row.names = FALSE)

driver_pyclone_df <-  driver_pyclone[,-6] %>%
  melt(id.vars = c("CHROM","POS","Hugo_Symbol"),
       measure.vars = c("MLS_1","RMLS_1"))

colnames(driver_pyclone_df)[4:5] <- c("Sample","CCF_value")


write.table(driver_pyclone_df,"pyclone_driver_df_p1.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#######################################################visualization of pyclone result 

library(ggrepel)

library(ggplot2)

library(ggpubr)

options(ggrepel.max.overlaps = Inf)

library(tibble)

scatter_for_py <- cluster_clonevol_df

scatter_for_py <- column_to_rownames(scatter_for_py,"mutation_id")

scatter_for_py$cluster <- as.factor(scatter_for_py$cluster)


table(scatter_for_py$cluster)


png(filename="pyclone_p1_driver_summary.png",width=1000,height=800,unit="px",bg="transparent",res = 120)

ggplot(scatter_for_py, aes(x =  MLS_1, y = RMLS_1)) +
  geom_jitter(aes(color = cluster), size = 2, alpha = 0.5, width = 0.01) +
  scale_color_discrete(guide = guide_legend(title = "Cluster"),
                       labels= c("1 (n=71)","2 (n=8)","3 (n=4)") ) + coord_equal() + theme_bw(base_size = 18) + theme(legend.position = "right")+geom_text_repel(data =scatter_for_py[driver_pyclone$Hugo_Symbol,],aes(label = driver_pyclone$Hugo_Symbol),
                                                                                 family = "Poppins",
                                                                                                                                                                 size = 3,
                                                                                                                                                                 min.segment.length = 0, 
                                                                                                                                                                 seed = 42, 
                                                                                                                                                                 box.padding = 0.5,
                                                                                                                                                                 max.overlaps = Inf,
                                                                                                                                                                 arrow = arrow(length = unit(0.010, "npc")),
                                                                                                                                                                 nudge_x = .15,
                                                                                                                                                                 nudge_y = .5,
                                                                                                                                                                 color = "grey20")+ expand_limits(x=c(0,100), y=c(0, 100))



dev.off()


pyclone_df <- cluster_clonevol_df %>%
  group_by(cluster) %>% 
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  melt(id.vars = c("mutation_id","cluster"),
       measure.vars = c("MLS_1","RMLS_1"))%>%
  mutate(cluster = as.factor(cluster)) %>%
  mutate(labels = ifelse(variable == "MLS_1",mutation_id, NA))

pyclone_driver_df <- subset(pyclone_df,pyclone_df$mutation_id %in% driver_pyclone$Hugo_Symbol)


png(filename="pyclone_p1_driver_point.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


 ggplot(pyclone_driver_df, aes(x=variable, y=value, group=mutation_id)) +
  geom_jitter(aes(colour=cluster), size=4.5, position=position_dodge(width=0.1)) + geom_line(alpha = 0.5, position=position_dodge(width=0.1)) +
   geom_text_repel(aes(label =labels, colour=cluster), direction = "y", nudge_x = -1,show.legend = FALSE,size = 2) + 
   xlab('Sample') +
  ylab('CCF value') +
  theme_bw()
 
 dev.off()
  
