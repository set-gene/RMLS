SE <- read.table("SE.MATS.JC.txt",sep = "\t", header = T)

SE <- SE[,c("geneSymbol","chr","strand","PValue","FDR","IncLevelDifference")]

cut_FDR <- 0.05

cut_lncd<- 0.05


SE_res<- SE[which(SE$FDR < cut_FDR & SE$IncLevelDifference > cut_lncd),]

SE_res$ASevent <- rep(c("SE"),len = nrow(SE_res))


##########################################################################
A3SS <- read.table("A3SS.MATS.JC.txt",sep = "\t", header = T)

A3SS <- A3SS[,c("geneSymbol","chr","strand","PValue","FDR","IncLevelDifference")]


A3SS_res<- A3SS[which(A3SS$FDR < cut_FDR & A3SS$IncLevelDifference > cut_lncd),]

A3SS_res$ASevent <- rep(c("A3SS"),len = nrow(A3SS_res))

#########################################################################

A5SS <- read.table("A5SS.MATS.JC.txt",sep = "\t", header = T)

A5SS <- A5SS[,c("geneSymbol","chr","strand","PValue","FDR","IncLevelDifference")]

A5SS_res<- A5SS[which(A5SS$FDR < cut_FDR & A5SS$IncLevelDifference > cut_lncd),]

A5SS_res$ASevent <- rep(c("A5SS"),len = nrow(A5SS_res))

##########################################################################

MXE <- read.table("MXE.MATS.JC.txt",sep = "\t", header = T)

MXE <- MXE[,c("geneSymbol","chr","strand","PValue","FDR","IncLevelDifference")]

MXE_res<- MXE[which(MXE$FDR < cut_FDR & MXE$IncLevelDifference > cut_lncd),]

MXE_res$ASevent <- rep(c("MXE"),len = nrow(MXE_res))

#########################################################################
RI <- read.table("RI.MATS.JC.txt",sep = "\t", header = T)

RI <- RI[,c("geneSymbol","chr","strand","PValue","FDR","IncLevelDifference")]

RI_res<- RI[which(RI$FDR < cut_FDR & RI$IncLevelDifference > cut_lncd),]

RI_res$ASevent <- rep(c("RI"),len = nrow(RI_res))

######################################################################

ASE_ALL <- rbind(A3SS_res,A5SS_res,MXE_res,RI_res,SE_res)

ASE_ALL <- ASE_ALL[,c(1,2,3,7,4,5)]

names(ASE_ALL)[1] <- c("Gene")

write.table(ASE_ALL, "ASE_RMLSFDR(0.05)t.txt", sep = "\t")

write.csv(x = ASE_ALL,"ASE_RMLSFDR(0.05).csv")
