#!/usr/bin/env bash
# Stream the 12.5 GB similarity.tsv.gz once and emit per-community within-community
# Jaccard distances for the 8 target communities (C1, C2, C3, C5, C6, C7, C11, C12).
# Output: paper_prep/figures/ed2/within_community_jaccard_<C>.tsv (cols: a, b, jaccard)
set -euo pipefail

ASSIGN=/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv
SIM=/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.similarity.tsv.gz
OUT=$(dirname "$0")

# Build seq -> arm-leiden community via 2-step lookup later; for ED2b we need
# the *arm-leiden* community labels (C1..C15), but the assignments file holds
# *seq-leiden* community ids in column 2. We need an arm -> arm-community map
# joined onto the seq table.
ARM=/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv

# Build per-sequence arm-leiden community label map
awk -F'\t' 'NR==FNR && FNR>1 {arm2c[$1]=$2; next} FNR>1 {if ($7 in arm2c) print $1"\t"arm2c[$7]}' \
  "$ARM" "$ASSIGN" > "$OUT/seq_to_arm_community.tsv"

echo "seq->arm_community map: $(wc -l < "$OUT/seq_to_arm_community.tsv") rows"

# Now stream similarity.tsv.gz, emit one row per pair where both seqs share
# an arm-community AND that community is in target list.
TARGETS="C1 C2 C3 C5 C6 C7 C11 C12"

zcat "$SIM" | \
awk -F'\t' -v map="$OUT/seq_to_arm_community.tsv" -v outdir="$OUT" -v targets="$TARGETS" '
BEGIN {
  while ((getline line < map) > 0) {
    n = split(line, f, "\t");
    if (n >= 2) seq2c[f[1]] = f[2];
  }
  close(map);
  ntgt = split(targets, T, " ");
  for (i = 1; i <= ntgt; i++) {
    keep[T[i]] = 1;
    fname[T[i]] = outdir "/within_community_jaccard_" T[i] ".tsv";
    print "a\tb\tjaccard" > fname[T[i]];
  }
}
NR == 1 { next }
{
  if (!($1 in seq2c) || !($2 in seq2c)) next;
  ca = seq2c[$1]; cb = seq2c[$2];
  if (ca != cb) next;
  if (!(ca in keep)) next;
  print $1 "\t" $2 "\t" $6 >> fname[ca];
}'

for c in $TARGETS; do
  f="$OUT/within_community_jaccard_${c}.tsv"
  if [[ -f "$f" ]]; then
    echo "$c: $(($(wc -l < "$f") - 1)) pairs"
  fi
done
