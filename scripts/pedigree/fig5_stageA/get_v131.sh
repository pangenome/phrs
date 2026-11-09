#!/bin/bash
# Download WashU pedigree v1.3.1 raw, bgzip, index. No header edits.
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies
U=https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.3.1
cd "$D"
for s in PAN010 PAN011 PAN027 PAN028; do
  for h in 1 2; do
    f=${s}_hap${h}.polished.v1.3.1.fasta
    out=${s}_hap${h}.v1.3.1.raw.fa.gz
    if [ -s "$out" ]; then echo "SKIP $out"; continue; fi
    echo "GET $f"
    curl -fsSL --retry 5 --retry-delay 10 -o "$f" "$U/$f"
    echo "BGZIP $out"
    $HOME/.local/bin/bgzip -@ 8 -c "$f" > "$out"
    rm -f "$f"
    echo "DONE $out"
  done
done
echo ALLDONE
