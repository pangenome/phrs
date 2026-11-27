#!/bin/bash
# Rest of calibrate_runs.sh after its run stopped at hg002_hic_mapq2 (an end without simulated cis contacts, now handled
# by calibrate_artefact.py): Hi-C MAPQ >= 2 (sim MAPQ >= 2 mcools already built by calibrate_runs.sh), Pore-C and CiFi.
set -euxo pipefail
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; O=$Q/calib
run () { $PY $Q/calibrate_artefact.py $Q/$1 $Q/$2 $3 $A/$4 $O/$5; }
run hg002_null.domains.tsv null_pk40/sim.mapq2.mcool $Q/mapq/HG002.mapq2.mcool hg002_phr_pair_correlation.tsv hg002_hic_mapq2
run chm13_null.domains.tsv chm13_null_pk40/sim.mapq2.mcool $Q/mapq/CHM13.mapq2.mcool chm13_phr_pair_correlation.tsv chm13_hic_mapq2
run hg002_null.domains.tsv porec_null_rk20/sim.mq1.mcool $HC/HG002/porec_output/porec.mcool hg002_phr_pair_correlation.tsv hg002_porec_mq1
run hg002_null.domains.tsv cifi_null_rk20/sim.mq1.mcool $HC/HG002/cifi_output/cifi.mcool hg002_phr_pair_correlation.tsv hg002_cifi_mq1
