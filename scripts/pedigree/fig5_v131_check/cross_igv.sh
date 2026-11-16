#!/bin/bash
# Extract 9q/3q reads from cross BAMs
set -uo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
O=check/ont/igv; mkdir -p $O
for v in v11 v131; do
  [ $v = v11 ] && { RS="PAN027#2#chr9:136168001-136188000"; RD="PAN011#2#chr9:133702042-133722041"; } || { RS="PAN027#2#chr9:136175732-136195731"; RD="PAN011#2#chr9:133702516-133722515"; }
  for b in ont/SRR378260*.PAN011.onPAN027.$v.bam; do $S view -F 0x904 $b "$RS" | cut -f1-6,10 || echo "ERR $b" >&2; done > $O/dad_on_self.$v.sam
  for b in ont/SRR378260*.PAN027.onPAN011.$v.bam; do $S view -F 0x904 $b "$RD" | cut -f1-6,10 || echo "ERR $b" >&2; done > $O/self_on_dad.$v.sam
  echo "$v dad_on_self $(wc -l < $O/dad_on_self.$v.sam) reads, BAMs: $(ls ont/SRR378260*.PAN011.onPAN027.$v.bam | xargs -n1 basename | cut -d. -f1 | tr '\n' ' ')"
  echo "$v self_on_dad $(wc -l < $O/self_on_dad.$v.sam) reads, BAMs: $(ls ont/SRR378260*.PAN027.onPAN011.$v.bam | xargs -n1 basename | cut -d. -f1 | tr '\n' ' ')"
done
ls $O/*.sam
