#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05
# Own reads: contig or assembly?
# args: WASHU, SAMPLE, VERSION, ASM_WINDOW (contig:start-end), CONTIG_FASTA, CONTIG_NAME, PICK_REGIONS (comma), OUT_TSV
# HiFi (hifi/bam) and ONT (ont_full) reads of SAMPLE from PICK_REGIONS, each aligned to (a) the assembly window and
# (b) the contig, competing in one reference. Per read: target, MAPQ, NM rate. Summary per set: reads to each, median NM,
# and reads close to it (NM rate <= 0.003 HiFi, <= 0.04 ONT: own-locus reads, not paralogs).
set -euo pipefail
D=$1; SM=$2; V=$3; AW=$4; CF=$5; CN=$6; PICK=$7; O=$8; cd $D
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; MM=$D/scripts_tmp/minimap2-2.31/minimap2
W=/scratch/$USER.cva.$SLURM_JOB_ID; mkdir -p $W; trap 'rm -rf $W' EXIT
log () { echo "[$(date -Is)] $*"; }
H=${AW%%#chr*}; $S faidx norm_$V/$H.fa.gz "$AW" | sed '1s/.*/>assembly/' > $W/ref.fa; $S faidx $CF "$CN" | sed '1s/.*/>contig/' >> $W/ref.fa
log "$SM $V: assembly $AW vs $CN; reads from $PICK"
printf "set\treads\tto_assembly\tto_contig\tnm_rate_med_assembly\tnm_rate_med_contig\tclose_assembly\tclose_contig\tclose_max\n" > $O
for set in hifi ont; do
  if [ $set = hifi ]; then B=$(ls hifi/bam/*.$SM.$V.bam); X=map-hifi; else B=$(ls ont_full/SRR*.$SM.$V.bam); X=map-ont; fi
  for b in $B; do $S view -F 0x900 $b $(echo $PICK | tr ',' ' '); done | awk -F'\t' '!s[$1]++ {print "@"$1"\n"$10"\n+\n"$11}' > $W/$set.fq
  $MM -ax $X --secondary=no -t 8 $W/ref.fa $W/$set.fq 2> /dev/null | $S view -F 0x904 -q 20 - \
    | awk -F'\t' -v s=$set -v OFS='\t' '{nm=0; for (i=12;i<=NF;i++) if ($i ~ /^NM:i:/) nm=substr($i,6); al=0; c=$6; while (match(c,/[0-9]+[M=XI]/)) {al+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; print s, $3, (al ? nm/al : 1)}' > $W/$set.tsv
  python3 - $W/$set.tsv $set $(awk 'NR%4==2' $W/$set.fq | wc -l) >> $O <<'PY'
import sys, statistics as st
f, s, n = sys.argv[1], sys.argv[2], sys.argv[3]
r = [l.split() for l in open(f)]
a = [float(x[2]) for x in r if x[1] == "assembly"]; c = [float(x[2]) for x in r if x[1] == "contig"]
m = 0.003 if s == "hifi" else 0.04  # close = own-locus reads (HiFi / ONT error level)
print(f"{s}\t{n}\t{len(a)}\t{len(c)}\t{st.median(a) if a else 'NA'}\t{st.median(c) if c else 'NA'}\t{sum(x <= m for x in a)}\t{sum(x <= m for x in c)}\t{m}")
PY
done
cat $O
