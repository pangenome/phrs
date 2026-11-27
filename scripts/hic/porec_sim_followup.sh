#!/bin/bash
# After a Pore-C null simulation (sim_porec_null.sbatch) in SIMDIR: the paper's contact statistics on the simulated
# matrices (parse2 --min-mapq 1 and 0) at 50 and 20 kb, the depth calibration of the paper's Pore-C matrix against the
# --min-mapq 1 simulation (calib/hg002_porec_mq1_TAG), and the paper's statistics on the paper's matrix minus the
# expected artefact (corrected_res/hg002_porec_TAG). The --min-mapq 0 comparison waits for the Pore-C rerun.
# usage: bash porec_sim_followup.sh SIMDIR TAG
set -euxo pipefail
D=$1; T=$2; Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; DOM=$Q/hg002_null.domains.tsv; FAI=$HC/HG002/HG002.fa.gz.fai
for mq in 1 0; do for r in 50000 20000; do bash $Q/contact_stats.sbatch --sample HG002_porec --matrix $D/sim.mq$mq.mcool --resolution $r --fai $FAI --out $D/stats_mq${mq}_$r; done; done
$PY $Q/calibrate_artefact.py $DOM $D/sim.mq1.mcool $HC/HG002/porec_output/porec.mcool $A/hg002_phr_pair_correlation.tsv $Q/calib/hg002_porec_mq1_$T
bash $Q/corrected_res_stats.sbatch --sample HG002_porec --domains $DOM --sim $D/sim.mq1.mcool --real $HC/HG002/porec_output/porec.mcool --fai $FAI --out $Q/corrected_res/hg002_porec_$T
