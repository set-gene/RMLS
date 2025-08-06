#!/bin/bash
cnvkit.py export seg *.cns -o MLS_gistic.seg
sed "s/.call.b//g" PT.gistic.seg > PT.gistic.sm.seg
sed '/chrM/d' PT.gistic.sm.seg > PT.gistic.chr.seg
