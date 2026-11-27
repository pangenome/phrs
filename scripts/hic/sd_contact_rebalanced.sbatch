#!/bin/bash
#SBATCH -c 4
#SBATCH -t 12:00:00
# R1-M5 on unique reads with PHR bins kept: the Hi-C MAPQ >= 2 matrices rebalanced without the low-coverage filter
# (cooler balance --force --mad-max 0, hicexplorer 3.7.4 cooler, as rebalance_stats.sbatch; the default filter masks
# 19 of 38 CHM13 PHR windows at 50 kb), then sd_contact_correlation.py as sd_contact_runs.sh, at 50 and 20 kb, for CHM13
# and HG002. The rebalanced matrices stay in node /scratch; the outputs go to matched_sd/contacts/<g>_hic_mapq2_madmax0_<r>.
# usage: sbatch sd_contact_rebalanced.sbatch
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
M=/moosefs/guarracino/HPRCv2/PHR_III/matched_sd; Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; SC=/moosefs/guarracino/HPRCv2/scripts/community
C=/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/cooler; PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3
mkdir -p /scratch/$USER/tmp/phrs; W=$(mktemp -d /scratch/$USER/tmp/phrs/sdreb.XXXX); trap 'rc=$?; rm -rf $W; echo "exit status $rc at line ${LINENO}"' EXIT
for g in chm13 hg002; do G=$(echo $g | tr a-z A-Z); sim=$(ls $M/pggb_${g}_cen/*.smooth.final.similarity.tsv.gz)
  for r in 50000 20000; do
    $C cp $Q/mapq/$G.mapq2.mcool::resolutions/$r $W/m.mcool::resolutions/$r
    $C balance --force --mad-max 0 $W/m.mcool::resolutions/$r
    $PY $M/sd_contact_correlation.py $SC $M/${g}_msd_cen_c0.5.pairs.tsv $M/${g}_msd_cen_c0.5.loci.pansn.fa.gz.names.tsv $sim $W/m.mcool $r $M/contacts/${g}_hic_mapq2_madmax0_$r
    rm -f $W/m.mcool
  done
done
echo "[$(date -Is)] done"
