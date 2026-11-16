#!/bin/bash
# Download WashU pedigree v1.1 raw, bgzip, index. No header edits.
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11
U=https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1
mkdir -p "$D"; cd "$D"
names="PAN010_hap1_HiFi_element_final_hap1 PAN010_hap2_HiFi_element_final_hap2 PAN011_hap1_HiFi_element_final_XY_hap1 PAN011_hap2_HiFi_element_final_XY_hap2 PAN027_mat_HiFi_element_final_mat PAN027_pat_HiFi_element_final_pat PAN028_hap1_HiFi_element_final_hap1 PAN028_hap2_HiFi_element_final_hap2"
for n in $names; do
  f=$n.polished.fasta
  out=$n.v1.1.raw.fa.gz
  if [ -s "$out" ]; then echo "SKIP $out"; continue; fi
  echo "GET $f"
  curl -fsSL --retry 5 --retry-delay 10 -o "$f" "$U/$f"
  $HOME/.local/bin/bgzip -@ 8 -c "$f" > "$out"
  rm -f "$f"
  echo "DONE $out"
done
echo ALLDONE_V11
