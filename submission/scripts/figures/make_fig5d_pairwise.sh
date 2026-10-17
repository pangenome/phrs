#!/usr/bin/env bash
# Fig5D - pairwise alignments of the sequences involved in the chr9q/chr3q case.
# Query: SELF.p.chr9q (PAN027 paternal chr9q, terminal 60 kb, telomere-first).
# Targets: DAD.h2.chr9q (PAN011 hap2 chr9q, transmitted to SELF),
#          DAD.h2.chr3q (PAN011 hap2 chr3q, the source),
#          DAD.h1.chr9q (negative control: his other chr9q has no module).
# Requires: the four 60 kb terminal windows + FastGA + pafplot.
set -euo pipefail
V11=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11
OUT=${1:-/tmp/fig5d}; mkdir -p "$OUT"
# 1. assemble the four telomere-first 60 kb windows
samtools faidx $V11/PAN027_pat_HiFi_element_final_pat.v1.1.raw.fa.gz PAN027.chr9.paternal:136149150-136209149 \
  | awk 'NR==1{print ">SELF.p.chr9q"} NR>1{print}' > $OUT/w.fa
revcomp () { awk 'NR==1{print;next}{print}' | python3 -c '
import sys
c=str.maketrans("ACGTNacgtn","TGCANtgcan")
l=sys.stdin.read().split("\n"); print(l[0])
s="".join(l[1:]).upper().translate(c)[::-1]
for i in range(0,len(s),80): print(s[i:i+80])'; }
for spec in "PAN011_hap2_HiFi_element_final_XY_hap2.v1.1.raw.fa.gz PAN011.chr9.haplotype2:133684562-133744561 DAD.h2.chr9q" \
            "PAN011_hap2_HiFi_element_final_XY_hap2.v1.1.raw.fa.gz PAN011.chr3.haplotype2:205518941-205578940 DAD.h2.chr3q" \
            "PAN011_hap1_HiFi_element_final_XY_hap1.v1.1.raw.fa.gz PAN011.chr9.haplotype1:136041079-136101078 DAD.h1.chr9q"; do
  set -- $spec
  samtools faidx $V11/$1 $2 | revcomp | sed "1s/.*/>$3/" >> $OUT/w.fa
done
# 2. sensitive, unfiltered alignments (short matches retained: -s200, high-freq seeds: -f8192)
FastGA -T8 -f8192 -s200 -c70 -i0.6 -paf $OUT/w.fa $OUT/w.fa > $OUT/pairwise.paf
# 3. plot: identity colouring on a 98-100% scale; tract boxed in each cell, no label text
pafplot $OUT/pairwise.paf --pdf --color-by-identity --color-scale 0.98:1.0 -s 2400 \
  --highlight "$(dirname "$0")/../../../fig/MainFigures/Fig5D_highlights.tsv" -o $OUT/Fig5D_pairwise_9q3q
echo "wrote $OUT/Fig5D_pairwise_9q3q.pdf"
