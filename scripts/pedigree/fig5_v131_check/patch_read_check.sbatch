#!/bin/bash
#SBATCH -p tux
#SBATCH -w tux05
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 12:00:00
# DAD reads on original and patched DAD
# args: WASHU, PATCH_DIR (OUT_DIR of patch_dad_chr9.sbatch), [VERSION] (default v11)
# Maps the DAD ONT reads of chr9q h2, chr1p h1, chr3q h1/h2 (own BAMs) and SELF.p chr9q (cross BAMs)
# to the original and the patched DAD diploid, then counts wrong positions and depth.
# out: PATCH_DIR/check/read_check.<v>.tsv
set -euo pipefail
D=$1; OUT=$2; V=${3:-v11}; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
M=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2
T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.patchchk.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
DC=PAN011#2#chr9; SC=PAN027#2#chr9
LD=$(grep -P "^$DC\t" norm_$V/PAN011#2.fa.gz.fai | cut -f2); LS=$(grep -P "^$SC\t" norm_$V/PAN027#2.fa.gz.fai | cut -f2)
P=$(awk -F'\t' 'NR == 2 {print $3}' $OUT/patch/patch.$V.tsv)
log "host $(hostname), $V, kept DAD $DC 1-$P (patch.$V.tsv)"
ls ont/SRR*.PAN011.onPAN027.$V.bam > $W/cross.txt
rm -f $OUT/check/read_check.$V.tsv; mkdir -p $OUT/check
for spec in "own:$DC:$((LD - 400000))-$LD" "own:PAN011#1#chr1:1-300000" "own:PAN011#1#chr3:END" "own:PAN011#2#chr3:END" "cross:$SC:$((LS - 400000))-$LS"; do
  k=${spec%%:*}; r=${spec#*:}; C=${r%%:*}
  if [ "${r#*:}" = END ]; then L=$(grep -P "^$C\t" norm_$V/${C%#chr*}.fa.gz.fai | cut -f2); r="$C:$((L - 300000))-$L"; fi
  if [ $k = own ]; then B=$(ls ont/SRR*.PAN011.$V.bam); else B=$(cat $W/cross.txt); fi
  for b in $B; do $S view -F 0x904 $b "$r"; done
done | awk -F'\t' '!s[$1]++ {print ">"$1"\n"$10}' > $W/reads.fa  # primary records keep the full read
log "reads.fa: $(grep -c '>' $W/reads.fa) reads, $(grep -v '>' $W/reads.fa | tr -d '\n' | wc -c) bp"
zcat norm_$V/PAN011#1.fa.gz norm_$V/PAN011#2.fa.gz > $W/orig.fa
zcat norm_$V/PAN011#1.fa.gz $OUT/norm_${V}_patched/PAN011#2.fa.gz > $W/patched.fa
pile () { python3 -c '
import re, sys
for l in sys.stdin:
    f = l.rstrip("\n").split("\t"); b = re.sub(r"\^.|\$", "", f[4]); i, s = 0, []
    while i < len(b):
        c = b[i]
        if c in "+-":
            m = re.match(r"[+-](\d+)", b[i:]); i += len(m.group(0)) + int(m.group(1)); continue
        s.append(c); i += 1
    nref = sum(c in ".," for c in s); alts = {x: sum(c.upper() == x for c in s) for x in "ACGT"}
    a = max(alts, key=alts.get); print(f"{f[1]}\t{len(s)}\t{nref}\t{a}\t{alts[a]}\t{f[2].upper()}")
'; }
for asm in orig patched; do
  $M -ax map-ont -I 8G -t $T $W/$asm.fa $W/reads.fa 2> $W/mm_$asm.log | $S sort -@ 8 -o $W/$asm.bam || { tail -5 $W/mm_$asm.log; exit 1; }
  log "$asm: $($S view -c -F 0x904 $W/$asm.bam) primary records"; $S index $W/$asm.bam; $S faidx $W/$asm.fa
  L=$(grep -P "^$DC\t" $W/$asm.fa.fai | cut -f2)
  for spec in "chr9q_h2_end:$DC:$((L - 330000))-$L" "chr1p_h1_start:PAN011#1#chr1:1-330000"; do
    n=${spec%%:*}; r=${spec#*:}
    $S mpileup -B -Q 0 -q 20 -d 100000 -f $W/$asm.fa -r "$r" $W/$asm.bam 2> /dev/null | pile > $W/pile.$asm.$n.tsv
    $S depth -a -r "$r" $W/$asm.bam > $W/depth0.$asm.$n.tsv
    python3 - $W/pile.$asm.$n.tsv $W/depth0.$asm.$n.tsv $asm $n $P $L <<'PY' | tee -a $OUT/check/read_check.$V.tsv
import sys, statistics as st
pile, dep, asm, n, P, L = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
# zones: chr9q end = after anchor vs 100-200 kb before anchor; chr1p = copy 8-61,926 vs 162-262 kb
z = (P + 1, L, P - 200000, P - 100000) if n.startswith("chr9") else (8, 61926, 162000, 262000)
rows = [l.rstrip("\n").split("\t") for l in open(pile)]
wrong = [(int(f[0]), f[5]) for f in rows if int(f[2]) + int(f[4]) >= 3 and int(f[4]) / (int(f[2]) + int(f[4])) >= 0.8]
open(f"{pile}.wrong", "w").writelines(f"{f[0]}\t{f[5]}\t{f[3]}\t{f[4]}\t{f[2]}\n" for f in rows if int(f[2]) + int(f[4]) >= 3 and int(f[4]) / (int(f[2]) + int(f[4])) >= 0.8)  # pos, asm base, read base, n read base, n asm base
d = {int(f[1]): int(f[2]) for f in (l.split("\t") for l in open(dep))}
inz = [d.get(i, 0) for i in range(z[0], z[1] + 1)]; ctl = [d.get(i, 0) for i in range(z[2], z[3] + 1)]
print(f"{asm}\t{n}\tzone={z[0]}-{z[1]}\twrong_positions_MAPQ20={sum(z[0] <= p <= z[1] and b != 'N' for p, b in wrong)}\tat_N_bases={sum(z[0] <= p <= z[1] and b == 'N' for p, b in wrong)}\twrong_positions_window={sum(b != 'N' for p, b in wrong)}\tdepth_all_zone_median={st.median(inz)}\tdepth_all_control_median={st.median(ctl)}")
PY
    cp $W/pile.$asm.$n.tsv.wrong $OUT/check/wrong_positions.$asm.$n.$V.tsv
  done
done
log done
