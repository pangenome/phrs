#!/bin/bash
# Corrected DAD v1.3.1, add chr7 fix
# run on octopus01: bash g1_build_dad_v131b.sh ; out $D/g1/corrected_v131b/norm/PAN011#{1,2}.fa.gz
# Start: g1/corrected_v131_pol (chr9 fixed). v1.3.1 chr7 h2 976,489-1,220,704 = v1.1 (244,216 bp exact, same
# coordinates), HiFi depth 2 in [1,076,488, 1,120,704) as v1.1. We put in the same polished piece as v1.1:
# g1/corrected_v11_pol PAN011#2#chr7 [1,076,488, 1,120,971) (44,483 bp, hifiasm --ul contig with DAD reads only).
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
[ $# -eq 0 ] || { echo "ERROR no arguments expected, got $*"; exit 2; }
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; I=g1/corrected_v131_pol; O=g1/corrected_v131b; P=g1/corrected_v11_pol
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
log () { echo "[$(date -Is)] $*"; }
[ -s $O/norm/PAN011#2.fa.gz.fai ] && { log "exists $O/norm/PAN011#2.fa.gz, skip"; exit 0; }
mkdir -p $O/norm /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/v131b.XXXX)
log "check flanks: v1.1 piece source vs v1.3.1 input, 100 kb each side"
for r in 976489-1076488 1120705-1220704; do
  $S faidx $I/norm/PAN011#2.fa.gz "PAN011#2#chr7:$r" | tail -n +2 | tr -d '\n' | md5sum > $W/a.md5
  $S faidx norm_v11/PAN011#2.fa.gz "PAN011#2#chr7:$r" | tail -n +2 | tr -d '\n' | md5sum > $W/b.md5
  cmp -s $W/a.md5 $W/b.md5 || { echo "ERROR flank $r differs between v1.3.1 and v1.1"; exit 1; }
  log "flank $r identical ($(cut -c1-12 $W/a.md5))"
done
$S faidx $P/norm/PAN011#2.fa.gz "PAN011#2#chr7:1-161478685" | sed '1s/.*/>v11pol_chr7/' > $W/piece.fa; $S faidx $W/piece.fa
python3 scripts_tmp/g1_apply_patch.py $I/norm/PAN011#2.fa.gz $W/out.fa 'PAN011#2#chr7' 1076488 1120704 $W/piece.fa v11pol_chr7 1076488 1120971 + 2> $O/chr7.log
cat $O/chr7.log
$BGZ -@ 4 -c $W/out.fa > $O/norm/PAN011#2.fa.gz && $S faidx $O/norm/PAN011#2.fa.gz
for f in PAN011#1.fa.gz PAN011#1.fa.gz.fai PAN011#1.fa.gz.gzi; do cp $I/norm/$f $O/norm/; done
log "lengths: $(grep -P 'chr7\t|chr9\t' $O/norm/PAN011#2.fa.gz.fai | cut -f1,2 | tr '\t\n' '  ')"
rm -rf "$W"
log done
