#!/bin/bash
# After a Hi-C null simulation (sim_hic_null.sbatch) of SAMPLE in SIMDIR: placement rates by MAPQ (sim_hic_rates.py),
# the paper's contact statistics at 50 and 20 kb, the MAPQ >= 1 and >= 2 simulated matrices (hic_mapq_filter.sbatch),
# the depth calibration against the real matrices (all reads, MAPQ >= 1, MAPQ >= 2; calib/<sample>_hic{,_mapq1,_mapq2}_TAG)
# and the paper's statistics on the real matrices minus the expected artefact at 50 and 20 kb (all reads and MAPQ >= 1;
# corrected_res/<sample>_hic{,_mapq1}_TAG).
# usage: bash hic_sim_followup.sh HG002|CHM13 SIMDIR TAG
set -euxo pipefail
SM=$1; D=$2; T=$3; s=$(echo $SM | tr A-Z a-z)
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; DOM=$Q/${s}_null.domains.tsv; FAI=$HC/$SM/$SM.fa.gz.fai
$PY $Q/sim_hic_rates.v2.py $D/sim.allValidPairs $D/eval.reads.tsv $D/rates_by_mapq.tsv
for r in 50000 20000; do bash $Q/contact_stats.sbatch --sample $SM --matrix $D/sim.mcool --resolution $r --fai $FAI --out $D/stats_$r; done
bash $Q/hic_mapq_filter.sbatch --pairs $D/sim.allValidPairs --sizes $HC/$SM/chrom.sizes --thresholds "1 2" --out $D/sim
$PY $Q/calibrate_artefact.py $DOM $D/sim.mcool $HC/$SM/hicpro_output/cool/$SM.mcool $A/${s}_phr_pair_correlation.tsv $Q/calib/${s}_hic_$T
for t in 1 2; do $PY $Q/calibrate_artefact.py $DOM $D/sim.mapq$t.mcool $Q/mapq/$SM.mapq$t.mcool $A/${s}_phr_pair_correlation.tsv $Q/calib/${s}_hic_mapq${t}_$T; done
bash $Q/corrected_res_stats.sbatch --sample $SM --domains $DOM --sim $D/sim.mcool --real $HC/$SM/hicpro_output/cool/$SM.mcool --fai $FAI --out $Q/corrected_res/${s}_hic_$T
bash $Q/corrected_res_stats.sbatch --sample $SM --domains $DOM --sim $D/sim.mapq1.mcool --real $Q/mapq/$SM.mapq1.mcool --fai $FAI --out $Q/corrected_res/${s}_hic_mapq1_$T
