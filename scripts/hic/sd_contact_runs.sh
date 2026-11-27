#!/bin/bash
# R1-M5: similarity vs contact for PHR loci and matched non-subtelomeric SD loci, per genome and contact matrix.
set -euxo pipefail
M=/moosefs/guarracino/HPRCv2/PHR_III/matched_sd; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; SC=/moosefs/guarracino/HPRCv2/scripts/community
mkdir -p $M/contacts
run () { g=$1; lab=$2; mc=$3; r=$4
  sim=$(ls $M/pggb_${g}_cen/*.smooth.final.similarity.tsv.gz)
  $PY $M/sd_contact_correlation.py $SC $M/${g}_msd_cen_c0.5.pairs.tsv $M/${g}_msd_cen_c0.5.loci.pansn.fa.gz.names.tsv $sim $mc $r $M/contacts/${g}_${lab}_${r}; }
for r in 50000 20000; do
  run chm13 hic $HC/CHM13/hicpro_output/cool/CHM13.mcool $r
  run hg002 hic $HC/HG002/hicpro_output/cool/HG002.mcool $r
  run hg002 porec $HC/HG002/porec_output/porec.mcool $r
  run hg002 cifi $HC/HG002/cifi_output/cifi.mcool $r
done
run chm13 hic_mapq2 $Q/mapq/CHM13.mapq2.mcool 50000
run hg002 hic_mapq2 $Q/mapq/HG002.mapq2.mcool 50000
