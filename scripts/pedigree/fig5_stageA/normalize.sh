#!/bin/bash
# Normalize per-hap WashU FASTAs to PanSN SAMPLE#H#chrN
set -euo pipefail
B=$HOME/.local/bin
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
SRC=$1; OUT=$2
mkdir -p "$OUT"
for f in "$SRC"/*.raw.fa.gz; do
  b=$(basename "$f"); samp=${b%%_*}; stem=${b%%.*}; tag=$(echo "$stem" | cut -d_ -f2)
  case $tag in hap1|mat) h=1;; hap2|pat) h=2;; *) echo "BAD TAG [$tag] from [$b]"; exit 1;; esac
  o="$OUT/${samp}#${h}.fa.gz"
  [ -s "$o.fai" ] && { echo "SKIP $o"; continue; }
  zcat "$f" | awk -v p="${samp}#${h}#" '
    /^>/ { keep=0
           if (match($1, /chr[0-9XY]+/)) { c=substr($1,RSTART,RLENGTH); keep=1; print ">" p c; next } }
    keep {print}
  ' | $B/bgzip -@ 8 > "$o"
  $S faidx "$o"
  echo "WROTE $o ($(wc -l < "$o.fai") records)"
done
echo NORMDONE
