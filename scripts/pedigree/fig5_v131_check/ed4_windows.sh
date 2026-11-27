#!/bin/bash
# Extended Data Fig. 4 inputs: the six telomere-first 60 kb ends (q ends reverse-complemented), FASTA on stdout.
# usage: ssh octopus01 bash -s < ed4_windows.sh > data/ed4_pairwise_windows.fa; then FastGA -T8 -f8192 -s200 -c70 -i0.6 -pafx w.fa w.fa > data/ed4_pairwise.paf
#   (FAtoGDB and GIXmake in PATH) and ed4_pairwise_dotplot.py data/ed4_pairwise_windows.fa data/ed4_pairwise.paf 21149 41149 OUT
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
fin () { awk -F'\t' -v s=$1 '$1 == "v11" && $2 == s {print $3}' scripts_tmp/g1_finals.tsv; }
SELF=$(fin PAN027); DAD=$(fin PAN011)
w () { local L; L=$(awk -v c="$2" '$1 == c {print $2}' "$1.fai"); $S faidx -i "$1" "$2:$((L - 59999))-$L" | sed "1s/.*/>$3/"; echo "$3 $1 $2 length $L" >&2; }
w $SELF/PAN027#2.fa.gz 'PAN027#2#chr9' SELF.p.chr9q
w norm_v11/PAN011#2.fa.gz 'PAN011#2#chr9' DAD.h2.chr9q.released
w $DAD/PAN011#2.fa.gz 'PAN011#2#chr9' DAD.h2.chr9q.corrected
w $DAD/PAN011#2.fa.gz 'PAN011#2#chr3' DAD.h2.chr3q
w $DAD/PAN011#1.fa.gz 'PAN011#1#chr3' DAD.h1.chr3q
w $DAD/PAN011#1.fa.gz 'PAN011#1#chr9' DAD.h1.chr9q
