#!/bin/bash
# Depth calibration (calibrate_artefact.py) for Hi-C restricted to MAPQ >= 1 on both reads: the MAPQ >= 1 versions of
# the two Hi-C null simulations (hic_mapq_filter.sbatch as a plain script) vs the real MAPQ >= 1 matrices.
set -euxo pipefail
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; O=$Q/calib; mkdir -p $O
for g in null_pk40:HG002 chm13_null_pk40:CHM13; do d=${g%%:*}; s=${g##*:}
  bash $Q/hic_mapq_filter.sbatch --pairs $Q/$d/sim.allValidPairs --sizes $HC/$s/chrom.sizes --thresholds 1 --out $Q/$d/sim
done
$PY $Q/calibrate_artefact.py $Q/hg002_null.domains.tsv $Q/null_pk40/sim.mapq1.mcool $Q/mapq/HG002.mapq1.mcool $A/hg002_phr_pair_correlation.tsv $O/hg002_hic_mapq1
$PY $Q/calibrate_artefact.py $Q/chm13_null.domains.tsv $Q/chm13_null_pk40/sim.mapq1.mcool $Q/mapq/CHM13.mapq1.mcool $A/chm13_phr_pair_correlation.tsv $O/chm13_hic_mapq1
