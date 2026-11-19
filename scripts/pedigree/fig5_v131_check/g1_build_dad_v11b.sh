#!/bin/bash
# Corrected DAD v1.1, chr9 from own contig
# run on octopus01: bash g1_build_dad_v11b.sh ; out $D/g1/corrected_v11b/norm/PAN011#{1,2}.fa.gz
# Start: g1/corrected_v11_pol PAN011#2 (chr7 fixed; chr9 up to 133,594,274 = original assembly).
# chr9 [133,594,274, end) <- DAD-only hifiasm contig h2tg000001l [144,573, end): the junction is the end of a
# 120,784 bp exact match contig = assembly (no SELF anchor, no SELF frame).
set -euxo pipefail
trap 'echo "exit at line ${LINENO}"' EXIT
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; O=g1/corrected_v11b; mkdir -p $O/norm
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
log () { echo "[$(date -Is)] $*"; }
L=$(awk '$1=="PAN011#2#chr9" {print $2}' g1/corrected_v11_pol/norm/PAN011#2.fa.gz.fai)
CL=$(awk '$1=="loc.bp.hap2.p_ctg.h2tg000001l" {print $2}' g1/local_asm/chr9_h2/contigs.fa.fai 2>/dev/null || true)
[ -n "$CL" ] || { $S faidx g1/local_asm/chr9_h2/contigs.fa; CL=$(awk '$1=="loc.bp.hap2.p_ctg.h2tg000001l" {print $2}' g1/local_asm/chr9_h2/contigs.fa.fai); }
log "chr9 length before $L, contig length $CL"
python3 scripts_tmp/g1_apply_patch.py g1/corrected_v11_pol/norm/PAN011#2.fa.gz $O/PAN011#2.fa 'PAN011#2#chr9' 133594274 $L \
  g1/local_asm/chr9_h2/contigs.fa loc.bp.hap2.p_ctg.h2tg000001l 144573 $CL + 2> $O/chr9.log
cat $O/chr9.log
$BGZ -@ 4 -c $O/PAN011#2.fa > $O/norm/PAN011#2.fa.gz && $S faidx $O/norm/PAN011#2.fa.gz
for f in PAN011#1.fa.gz PAN011#1.fa.gz.fai PAN011#1.fa.gz.gzi PAN027#1.fa.gz PAN027#1.fa.gz.fai PAN027#1.fa.gz.gzi PAN027#2.fa.gz PAN027#2.fa.gz.fai PAN027#2.fa.gz.gzi; do cp g1/corrected_v11_pol/norm/$f $O/norm/; done
grep -P "chr7\t|chr9\t" $O/norm/PAN011#2.fa.gz.fai | cut -f1,2
log done
