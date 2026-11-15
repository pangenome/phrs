#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --exclude=octopus05
# Local HiFi reassembly of an error
# args: WASHU, SAMPLE, VERSION, ERROR (contig:start-end, 1-based), FLANK_BP, OUT_DIR
# 1. Window = ERROR +- FLANK on the faulty haplotype; its homologous window on the other haplotype from minimap2 asm5.
# 2. HiFi reads of SAMPLE (hifi/bam, primary) in both windows -> hifiasm 0.25.0 (system gcc build) contigs.
# 3. Contigs vs the faulty window and the other window (minimap2 asm5): span, identity; a contig "bridges" when it aligns to
#    both flanks of the faulty window (>= 20 kb each side, identity >= 0.999) with the error inside.
# out: OUT_DIR/contigs.fa, OUT_DIR/contigs_vs_windows.paf, OUT_DIR/summary.tsv
set -euo pipefail
D=$1; SM=$2; V=$3; ERR=$4; F=$5; O=$6; cd $D; mkdir -p $O
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
HA=$D/scripts_tmp/hifiasm-0.25.0/hifiasm; T=$SLURM_CPUS_PER_TASK
W=/scratch/$USER.hlasm.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
C=${ERR%%:*}; R=${ERR#*:}; A=${R%-*}; Z=${R#*-}; H=$(echo $C | cut -d'#' -f2); OH=$((3 - H)); OC=${C/\#$H\#/\#$OH\#}
L=$(awk -v c="$C" '$1==c {print $2}' norm_$V/$SM#$H.fa.gz.fai); WS=$((A > F ? A - F : 1)); WE=$((Z + F > L ? L : Z + F))
log "$SM $V error $ERR, window $C:$WS-$WE (length $L), other haplotype $OC, host $(hostname), hifiasm $($HA --version)"
$S faidx norm_$V/$SM#$H.fa.gz "$C:$WS-$WE" | sed '1s/.*/>faulty/' > $W/faulty.fa
$S faidx norm_$V/$SM#$OH.fa.gz "$OC" > $W/other_chr.fa
$MM -cx asm5 -t $T $W/other_chr.fa $W/faulty.fa 2> /dev/null | sort -k10,10nr > $W/proj.paf
read OS OE < <(awk -F'\t' '$12 >= 5 {s = (s == "" || $8 < s) ? $8 : s; e = ($9 > e) ? $9 : e} END {print s, e}' $W/proj.paf)
[ -n "${OS:-}" ] || { log "ERROR no homologous window on $OC"; exit 1; }
log "homologous window $OC:$OS-$OE ($(wc -l < $W/proj.paf) chains)"
$S faidx $W/other_chr.fa "$OC:$((OS + 1))-$OE" | sed '1s/.*/>other/' > $W/other.fa
for b in hifi/bam/*.$SM.$V.bam; do $S view -F 0x900 $b "$C:$WS-$WE" "$OC:$((OS + 1))-$OE"; done \
  | awk -F'\t' '!s[$1]++ {print "@"$1"\n"$10"\n+\n"$11}' > $W/reads.fq
log "reads $(awk 'NR%4==2' $W/reads.fq | wc -l), bp $(awk 'NR%4==2 {s+=length($0)} END {print s}' $W/reads.fq)"
cd $W; $HA -t $T -o loc reads.fq > hifiasm.log 2>&1 || { tail -5 hifiasm.log; exit 1; }
for g in loc.bp.hap1.p_ctg loc.bp.hap2.p_ctg loc.bp.p_ctg; do
  [ -s $g.gfa ] && awk -v p=$g '$1=="S" {print ">"p"."$2"\n"$3}' $g.gfa
done > $O/contigs.fa
log "contigs $(grep -c '>' $O/contigs.fa), bp $(grep -v '>' $O/contigs.fa | tr -d '\n' | wc -c)"
cat $W/faulty.fa $W/other.fa > $W/wins.fa
$MM -cx asm5 --eqx -t $T $W/wins.fa $O/contigs.fa 2> /dev/null > $O/contigs_vs_windows.paf
python3 - $O/contigs_vs_windows.paf $((A - WS)) $((Z - WS + 1)) $O/summary.tsv <<'PY'
import sys, collections
paf, ea, ez, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
r = collections.defaultdict(list)
for l in open(paf):
    f = l.rstrip("\n").split("\t"); r[f[0]].append((f[5], int(f[7]), int(f[8]), int(f[9]) / int(f[10]), int(f[1]), int(f[2]), int(f[3])))
with open(out, "w") as o:
    o.write("contig\tlength\ttarget\tt_start\tt_end\tidentity\tbridges_faulty\n")
    for c, hits in r.items():
        fl = [h for h in hits if h[0] == "faulty" and h[3] >= 0.999]
        left = any(h[1] <= ea - 20000 and h[2] >= ea for h in fl) or any(h[2] <= ea and ea - h[1] >= 20000 for h in fl)
        right = any(h[1] <= ez and h[2] >= ez + 20000 for h in fl) or any(h[1] >= ez and h[2] - h[1] >= 20000 for h in fl)
        for h in sorted(hits, key=lambda h: -(h[2] - h[1])):
            o.write(f"{c}\t{h[4]}\t{h[0]}\t{h[1]}\t{h[2]}\t{h[3]:.5f}\t{'yes' if left and right else 'no'}\n")
print(open(out).read())
PY
log done
