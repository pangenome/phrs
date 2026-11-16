#!/bin/bash
# Build Fig5 joint query/target FASTAs (mirrors prepare_inputs.py)
set -euo pipefail
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools
NORM=$1; OUT=$2
mkdir -p "$OUT"
# comparison_id  query_sample query_hap  target_sample
while read -r cid qs qh ts; do
  q="$OUT/${cid}.query.fa"; t="$OUT/${cid}.target.fa"
  if [ ! -s "$q.fai" ]; then
    zcat "$NORM/${qs}#${qh}.fa.gz" > "$q"; $S faidx "$q"
  fi
  if [ ! -s "$t.fai" ]; then
    # concat both parent haps, rewrite SAMPLE#H#name -> SAMPLE#joint#hH_name
    for h in 1 2; do
      zcat "$NORM/${ts}#${h}.fa.gz" | sed "s|^>${ts}#${h}#|>${ts}#joint#h${h}_|"
    done > "$t"
    $S faidx "$t"
  fi
  echo "BUILT $cid  query=$(wc -l < "$q.fai")  target=$(wc -l < "$t.fai")"
done <<TSV
PAN027pat_vs_PAN011_joint PAN027 2 PAN011
PAN027mat_vs_PAN010_joint PAN027 1 PAN010
PAN028mat_vs_PAN027_joint PAN028 1 PAN027
TSV
echo INPUTSDONE
