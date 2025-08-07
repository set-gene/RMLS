#!/bin/bash
oncodrivefml -i ./MLS_oncodrive_input.txt -e /SSD/jeeh99/tools/oncodrivefml2/resource/cds_hg38_chr.tsv -c /SSD/jeeh99/tools/oncodrivefml2/resource/oncodrivefml_v2.conf --sequencing wes -o ./MLS_re --cores 8

oncodrivefml -i ./RMLS_oncodrive_input.txt -e /SSD/jeeh99/tools/oncodrivefml2/resource/cds_hg38_chr.tsv -c /SSD/jeeh99/tools/oncodrivefml2/resource/oncodrivefml_v2.conf --sequencing wes -o ./RMLS_re --cores 8
