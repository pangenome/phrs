#!/bin/bash
# After porec_mapq_rerun.sbatch for ASSAY (porec or cifi): (1) the default branch (mq1) must equal the paper's matrix
# pixel by pixel at 5 kb (compare_cool.py); (2) the paper's statistics on mq1 and mq0 at 50 and 20 kb; (3) depth
# calibration of the mq0 matrix against the mq0 null simulation; (4) the statistics on the mq0 matrix minus the expected
# artefact (corrected_stats.sbatch as a plain script).
# Step 4 subtracts at 5 kb (corrected_stats.sbatch), which removes only part of the artefact; the jobs queued on
# 2026-11-26 ran it so, and corrected_res_stats.sbatch (50 and 20 kb) ran after them on the mq0 matrices.
# usage: bash rerun_followup.sh porec|cifi
set -euxo pipefail
AS=$1; Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; R=$Q/rerun_mapq/$AS; FAI=$HC/HG002/HG002.fa.gz.fai; SM=HG002_$AS
$PY $Q/compare_cool.py $R/mq1/$AS.mcool $HC/HG002/${AS}_output/$AS.mcool
for mq in 1 0; do for r in 50000 20000; do bash $Q/contact_stats.sbatch --sample $SM --matrix $R/mq$mq/$AS.mcool --resolution $r --fai $FAI --out $R/stats_mq${mq}_$r; done; done
$PY $Q/calibrate_artefact.py $Q/hg002_null.domains.tsv $Q/${AS}_null_rk20/sim.mq0.mcool $R/mq0/$AS.mcool $A/hg002_phr_pair_correlation.tsv $Q/calib/hg002_${AS}_mq0
bash $Q/corrected_stats.sbatch --sample $SM --domains $Q/hg002_null.domains.tsv --sim $Q/${AS}_null_rk20/sim.mq0.mcool --real $R/mq0/$AS.mcool --fai $FAI --out $Q/corrected/hg002_${AS}_mq0
