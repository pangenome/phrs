#!/bin/bash
# Rest of hic_sim_followup.sh HG002 null_wprof_pk40 wprof after it stopped at the MAPQ >= 2 calibration (pairs of two
# ends without a simulated scale, now bounded by the smallest and largest scale in calibrate_artefact.py).
set -euxo pipefail
SM=HG002; s=hg002; D=/moosefs/guarracino/HPRCv2/PHR_III/simhic/null_wprof_pk40; T=wprof
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3; DOM=$Q/${s}_null.domains.tsv; FAI=$HC/$SM/$SM.fa.gz.fai
$PY $Q/calibrate_artefact.py $DOM $D/sim.mapq2.mcool $Q/mapq/$SM.mapq2.mcool $A/${s}_phr_pair_correlation.tsv $Q/calib/${s}_hic_mapq2_$T
bash $Q/corrected_res_stats.sbatch --sample $SM --domains $DOM --sim $D/sim.mcool --real $HC/$SM/hicpro_output/cool/$SM.mcool --fai $FAI --out $Q/corrected_res/${s}_hic_$T
bash $Q/corrected_res_stats.sbatch --sample $SM --domains $DOM --sim $D/sim.mapq1.mcool --real $Q/mapq/$SM.mapq1.mcool --fai $FAI --out $Q/corrected_res/${s}_hic_mapq1_$T
