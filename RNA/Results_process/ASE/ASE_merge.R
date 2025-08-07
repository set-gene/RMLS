library(ggVennDiagram)

ASE_MR = read.csv("ASE_ALLtop(0.05)_MR.csv")

ASE_MR = ASE_MR[,-1]

ASE_NM = read.csv("ASE_ALLtop(0.05)_NM.csv")

ASE_NM = ASE_NM[,-1]

ASE_NR = read.csv("ASE_ALLtop(0.05)_NR.csv")

ASE_NR = ASE_NR[,-1]

ASE_MR_G = ASE_MR[,2]

ASE_MR_G = unique(ASE_MR_G)

ASE_NR_G = ASE_NR[,2]

ASE_NR_G = unique(ASE_NR_G)

ASE_NM_G = ASE_NM[,2]

ASE_NM_G = unique(ASE_NM_G)

dat_ASE = list(
  A = ASE_MR_G,
  B = ASE_NR_G,
  C = ASE_NM_G
)

library(ggplot2)

ggVennDiagram(dat_ASE,label_alpha = 0,
              category.names = c("MLS vs RMLS","NL vs RMLS","NL vs MLS"))+
                scale_fill_gradient(low = "blue", high = "red")+
                ggtitle("Comparison of ASE gene sets") 



library(ggVennDiagram)

library(ggpubr)

library(ggplot2)

library(officer)
library(rvg)

pvenn <- ggVennDiagram(dat_ASE, 
                       label_alpha = 0,
                       category.names = c("MLS vs RMLS","NL vs RMLS","NL vs MLS"),
                       set_color = c("#D55E00","#CC79A7","#1B9E77"),) +
  scale_x_continuous(expand = expansion(mult = .6)) + scale_fill_distiller(palette = "Reds", direction = 1)

#####################################
read_pptx() %>% 
  print(target = "./example.pptx")

read_pptx( "./example.pptx")


editable_graph <- rvg::dml(ggobj = pvenn)

read_pptx( "./example.pptx") %>%
  add_slide("Title Slide","Office Theme") %>%
  ph_with(editable_graph,
          location = ph_location(left=0, top=0,width = 6, height = 6, bg="transparent")) %>%
  print(target = "./example.pptx")
#####################################


library(ggvenn)

names(dat_ASE) <- c("MLS vs RMLS","NL vs RMLS","NL vs MLS")


png(filename="final_result/ASE_venn.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


ggvenn(dat_ASE, columns = c("MLS vs RMLS","NL vs RMLS","NL vs MLS"),
       stroke_size = 0.5,fill_color = c("blue","red","purple"),set_name_size = 5)+
  ggtitle("Comparison of ASE gene sets") + theme(plot.title = element_text(size = 15, face = "bold"))

dev.off()

ASE_MR$IncLevelDifference

cut_PSI <- 0.1

cut_fdr <- 0.05

ASE_MR_RMLS <- ASE_MR[which(ASE_MR$FDR < cut_fdr & (ASE_MR$IncLevelDifference<(-cut_PSI))),]

ASE_MR_RMLS_MICAL3 <- subset(ASE_MR_RMLS,ASE_MR_RMLS$geneSymbol == "MICAL3")

# Adjusted P values



png(filename="final_result/sig_ASE_MR.png",width=800,height=600,unit="px",bg="transparent",res = 150)



with(ASE_MR_RMLS , plot(IncLevelDifference, -log10(FDR), pch=20, main="Scatter plot for significant ASE", col='grey', xlab=bquote(~delta~PSI), ylab=bquote(~-log[10]~FDR),cex = 2))


intenteded_gene = c("MICAL3")


with(ASE_MR_RMLS_MICAL3, {
  points(IncLevelDifference, -log10(FDR), col="red", cex=2, lwd=2)
  text(IncLevelDifference, -log10(FDR), intenteded_gene, pos=2, col="red")
})
dev.off()
####################################################################
Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

intersect_genes <- data.frame(Intersect(dat_ASE))

library(dplyr)

ASE_MR_CD37_A3SS = read.csv("MLS_RMLS_CD37_A3SS.csv")

ASE_NM_CD37_A3SS = read.csv("NL_MLS_CD37_A3SS.csv")

ASE_CD37_A3SS = left_join(ASE_MR_CD37_A3SS,ASE_NM_CD37_A3SS, by = "event",keep = FALSE)

ASE_CD37_A3SS = ASE_CD37_A3SS[,-c(20:25)]

colnames(ASE_CD37_A3SS)[2:7] = c("MLS_1","MLS_2","MLS_3","MLS_4","MLS_5","MLS_6")

MLS_CD37_A3SS = ASE_CD37_A3SS[,2:7]

MLS_CD37_A3SS = data.frame(t(MLS_CD37_A3SS))

colnames(MLS_CD37_A3SS) = 'value'

MLS_CD37_A3SS$group = rep('MLS')

RMLS_CD37_A3SS = ASE_CD37_A3SS[,8:13]

RMLS_CD37_A3SS = data.frame(t(RMLS_CD37_A3SS))

colnames(RMLS_CD37_A3SS) = 'value'

RMLS_CD37_A3SS$group = rep('RMLS')

NL_CD37_A3SS = ASE_CD37_A3SS[,14:19]

NL_CD37_A3SS = data.frame(t(NL_CD37_A3SS))

colnames(NL_CD37_A3SS) = 'value'

NL_CD37_A3SS$group = rep('NL')

for_BAR <- rbind(MLS_CD37_A3SS, RMLS_CD37_A3SS)

for_BAR <- rbind(for_BAR, NL_CD37_A3SS)

colnames(for_BAR) <-c("Value", "Condition")

library(ggplot2)

library(ggsignif)

library(ggpubr)


for_BAR$Condition <- factor(for_BAR$Condition,levels = c("NL","MLS","RMLS"))


comparisons = list(c("MLS","RMLS"), c("RMLS","NL"), c("MLS","NL"))


png(filename="final_result/A3SS_PSI_value_CD37.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

ggboxplot(for_BAR, x="Condition", y="Value",color = "Condition",palette = c("green","blue", "red"),ylab = "A3SS PSI value of CD37") + stat_compare_means(method = "anova")

dev.off()


png(filename="final_result/A3SS_PSI_value_CD37_final.png",width=800,height=600,unit="px",bg="transparent",res = 150)

ggplot(for_BAR, aes(x=Condition, y=Value, fill=Condition)) + 
  geom_boxplot()+
  scale_y_continuous(name = "A3SS PSI value of CD37") +
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
  scale_fill_manual(values = c("#00B050","#1E1EB4","#DF1813"))+ 
  stat_compare_means(method = "anova")


dev.off()
#########################################################################################################

ASE_MR_CD37_RI = read.csv("MLS_RMLS_CD37_RI .csv")

ASE_NM_CD37_RI = read.csv("NL_MLS_CD37_RI .csv")

ASE_CD37_RI = left_join(ASE_MR_CD37_RI,ASE_NM_CD37_RI, by = "event",keep = FALSE)

ASE_CD37_RI= ASE_CD37_RI[,-c(20:25)]

colnames(ASE_CD37_RI)[2:7] = c("MLS_1","MLS_2","MLS_3","MLS_4","MLS_5","MLS_6")

MLS_CD37_RI = ASE_CD37_RI[,2:7]

MLS_CD37_RI = data.frame(t(MLS_CD37_RI))

colnames(MLS_CD37_RI) = 'value'

MLS_CD37_RI$group = rep('MLS')

RMLS_CD37_RI = ASE_CD37_RI[,8:13]

RMLS_CD37_RI = data.frame(t(RMLS_CD37_RI))

colnames(RMLS_CD37_RI) = 'value'

RMLS_CD37_RI$group = rep('RMLS')

NL_CD37_RI = ASE_CD37_RI[,14:19]

NL_CD37_RI = data.frame(t(NL_CD37_RI))

colnames(NL_CD37_RI) = 'value'

NL_CD37_RI$group = rep('NL')

for_BAR_RI <- rbind(MLS_CD37_RI, RMLS_CD37_RI)

for_BAR_RI <- rbind(for_BAR_RI, NL_CD37_RI)

colnames(for_BAR_RI) <-c("Value", "Condition")
library(ggplot2)

library(ggsignif)

library(ggpubr)


for_BAR_RI$Condition<- factor(for_BAR_RI$Condition,levels = c("NL","MLS","RMLS"))

comparisons = list(c("MLS","RMLS"), c("RMLS","NL"), c("MLS","NL"))


png(filename="final_result/RI_PSI_value_CD37.png",width=1000,height=800,unit="px",bg="transparent",res = 150)

ggboxplot(for_BAR_RI, x="Condition", y="Value",color = "Condition",palette = c("green","blue", "red"),ylab = "RI PSI value of CD37") + stat_compare_means(method = "anova")

dev.off()


png(filename="final_result/RI_PSI_value_CD37_final.png",width=800,height=600,unit="px",bg="transparent",res = 150)

ggplot(for_BAR_RI, aes(x=Condition, y=Value, fill=Condition)) + 
  geom_boxplot()+
  scale_y_continuous(name = "RI PSI value of CD37") +
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
  scale_fill_manual(values = c("#00B050","#1E1EB4","#DF1813"))+ 
  stat_compare_means(method = "anova")


dev.off()
