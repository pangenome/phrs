#!/bin/bash
# After the CHM13 Hi-C null simulation with the real read length (150 bp, chm13_null_rl150): placement rates by MAPQ
# (sim_hic_rates.py), the paper's contact statistics at 50 and 20 kb, the MAPQ >= 1 and >= 2 versions of the simulated
# matrix, and the depth calibration against the real CHM13 matrices (all reads, MAPQ >= 1, MAPQ >= 2).
set -euxo pipefail
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; D=$Q/chm13_null_rl150
$PY $Q/sim_hic_rates.v2.py $D/sim.allValidPairs $D/eval.reads.tsv $D/rates_by_mapq.tsv
for r in 50000 20000; do bash $Q/contact_stats.sbatch --sample CHM13 --matrix $D/sim.mcool --resolution $r --fai $HC/CHM13/CHM13.fa.gz.fai --out $D/stats_$r; done
bash $Q/hic_mapq_filter.sbatch --pairs $D/sim.allValidPairs --sizes $HC/CHM13/chrom.sizes --thresholds "1 2" --out $D/sim
$PY $Q/calibrate_artefact.py $Q/chm13_null.domains.tsv $D/sim.mcool $HC/CHM13/hicpro_output/cool/CHM13.mcool $A/chm13_phr_pair_correlation.tsv $Q/calib/chm13_hic_rl150
for t in 1 2; do $PY $Q/calibrate_artefact.py $Q/chm13_null.domains.tsv $D/sim.mapq$t.mcool $Q/mapq/CHM13.mapq$t.mcool $A/chm13_phr_pair_correlation.tsv $Q/calib/chm13_hic_mapq${t}_rl150; done
