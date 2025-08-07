library(maser)

library(rtracklayer)

path <- system.file("extdata", file.path("ASE_H_NL"),
                    package = "maser")

ris <- maser(path, c("NL", "RIS"), ftype = "JC")

ris

head(summary(ris, type = "SE")[, 1:8])

mls_filt <- filterByCoverage(ris, avg_reads = 5)

mls_filt

mls_top <- topEvents(mls_filt, fdr = 0.05, deltaPSI = 0.05)

mls_top@SE_stats

#############################################################################PARSINGSTAT

all_se <- as.data.frame(mls_top@SE_stats)

all_se$ASevent <- rep(c("SE"),len = nrow(all_se))

all_a3ss <- as.data.frame(mls_top@A3SS_stats)

all_a3ss$ASevent <- rep(c("A3SS"),len = nrow(all_a3ss))

all_a5ss <- as.data.frame(mls_top@A5SS_stats)

all_a5ss$ASevent <- rep(c("A5SS"),len = nrow(all_a5ss))

all_mxe <- as.data.frame(mls_top@MXE_stats)

all_mxe$ASevent <- rep(c("MXE"),len = nrow(all_mxe))

all_ri <- as.data.frame(mls_top@RI_stats)

all_ri$ASevent <- rep(c("RI"),len = nrow(all_ri))
##############################################################################
library(dplyr)

se_gene <- as.data.frame(mls_top@SE_events)

all_se <- left_join(all_se,se_gene,by="ID")

a3ss_gene <- as.data.frame(mls_top@A3SS_events)

all_a3ss <- left_join(all_a3ss,a3ss_gene,by="ID")

a5ss_gene <- as.data.frame(mls_top@A5SS_events)

all_a5ss <- left_join(all_a5ss,a5ss_gene,by="ID")

mxe_gene <- as.data.frame(mls_top@MXE_events)

all_mxe <- left_join(all_mxe,mxe_gene,by="ID")

ri_gene <- as.data.frame(mls_top@RI_events)

all_ri <- left_join(all_ri,ri_gene,by="ID")


###########################################################################


ASE_ALL <- rbind(all_a3ss,all_a5ss,all_mxe,all_ri,all_se)

ASE_ALL <- ASE_ALL[,c(6,7,2,3,4,5)]

ASE_High <- ASE_ALL[!(ASE_ALL$IncLevelDifference > 0 ),]

ASE_Low <- ASE_ALL[!(ASE_ALL$IncLevelDifference < 0 ),]

write.csv(x = ASE_ALL,"ASE_ALLtop(0.05).csv")

write.csv(x = ASE_RMLS,"ASE_RMLStop(0.05).csv")

write.csv(x = ASE_MLS,"ASE_MLStop(0.05).csv")
