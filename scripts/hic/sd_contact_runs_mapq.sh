#!/bin/bash
# R1-M5 follow-up: similarity vs contact for PHR loci and matched SD loci (sd_contact_correlation.py, as sd_contact_runs.sh)
# on the Hi-C matrices restricted to MAPQ >= 1 (50 and 20 kb) and MAPQ >= 2 (20 kb; 50 kb in sd_contact_runs.sh), and on
# the CiFi rerun with multi-mapping reads kept (parse2 --min-mapq 0, 50 and 20 kb).
set -euxo pipefail
M=/moosefs/guarracino/HPRCv2/PHR_III/matched_sd; Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; SC=/moosefs/guarracino/HPRCv2/scripts/community
run () { g=$1; lab=$2; mc=$3; r=$4
  sim=$(ls $M/pggb_${g}_cen/*.smooth.final.similarity.tsv.gz)
  $PY $M/sd_contact_correlation.py $SC $M/${g}_msd_cen_c0.5.pairs.tsv $M/${g}_msd_cen_c0.5.loci.pansn.fa.gz.names.tsv $sim $mc $r $M/contacts/${g}_${lab}_${r}; }
for r in 50000 20000; do
  run chm13 hic_mapq1 $Q/mapq/CHM13.mapq1.mcool $r
  run hg002 hic_mapq1 $Q/mapq/HG002.mapq1.mcool $r
  run hg002 cifi_mq0 $Q/rerun_mapq/cifi/mq0/cifi.mcool $r
done
run chm13 hic_mapq2 $Q/mapq/CHM13.mapq2.mcool 20000
run hg002 hic_mapq2 $Q/mapq/HG002.mapq2.mcool 20000
