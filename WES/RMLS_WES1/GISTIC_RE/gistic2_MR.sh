#!/bin/bash

gistic2 -b /HDD8T3/jeeh99/RMLS_WES1/GISTIC_RE/CNVkit_gistic_MLS/ -seg /HDD8T3/jeeh99/RMLS_WES1/CNVkit2/MLS_CNS/MLS_gistic.seg -refgene /SSD/jeeh99/gistic2/hg38.UCSC.add_miR.160920.refgene.mat -conf 0.99 >> ./MLS_gistic2.log 2>&1

gistic2 -b /HDD8T3/jeeh99/RMLS_WES1/GISTIC_RE/CNVkit_gistic_RMLS/ -seg /HDD8T3/jeeh99/RMLS_WES1/CNVkit2/RMLS_CNS/RMLS_gistic.seg -refgene /SSD/jeeh99/gistic2/hg38.UCSC.add_miR.160920.refgene.mat -conf 0.99 >> ./RMLS_gistic2.log 2>&1
