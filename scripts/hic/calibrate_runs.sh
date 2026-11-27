#!/bin/bash
# Depth calibration of the mapping artefact (calibrate_artefact.py) for each assay: null simulation vs real matrix.
# First the MAPQ >= 2 versions of the two Hi-C null simulations (hic_mapq_filter.sbatch as a plain script, as for the
# real MAPQ >= 2 matrices). Pore-C / CiFi: parse2 default (--min-mapq 1) simulation vs the paper's matrices.
set -euxo pipefail
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; O=$Q/calib; mkdir -p $O
for g in null_pk40:HG002 chm13_null_pk40:CHM13; do d=${g%%:*}; s=${g##*:}
  bash $Q/hic_mapq_filter.sbatch --pairs $Q/$d/sim.allValidPairs --sizes $HC/$s/chrom.sizes --thresholds 2 --out $Q/$d/sim
done
run () { $PY $Q/calibrate_artefact.py $Q/$1 $Q/$2 $3 $A/$4 $O/$5; }
run hg002_null.domains.tsv null_pk40/sim.mcool $HC/HG002/hicpro_output/cool/HG002.mcool hg002_phr_pair_correlation.tsv hg002_hic
run chm13_null.domains.tsv chm13_null_pk40/sim.mcool $HC/CHM13/hicpro_output/cool/CHM13.mcool chm13_phr_pair_correlation.tsv chm13_hic
run hg002_null.domains.tsv null_pk40/sim.mapq2.mcool $Q/mapq/HG002.mapq2.mcool hg002_phr_pair_correlation.tsv hg002_hic_mapq2
run chm13_null.domains.tsv chm13_null_pk40/sim.mapq2.mcool $Q/mapq/CHM13.mapq2.mcool chm13_phr_pair_correlation.tsv chm13_hic_mapq2
run hg002_null.domains.tsv porec_null_rk20/sim.mq1.mcool $HC/HG002/porec_output/porec.mcool hg002_phr_pair_correlation.tsv hg002_porec_mq1
run hg002_null.domains.tsv cifi_null_rk20/sim.mq1.mcool $HC/HG002/cifi_output/cifi.mcool hg002_phr_pair_correlation.tsv hg002_cifi_mq1
