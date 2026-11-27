#!/bin/bash
# Run contact_stats.sbatch (as a plain script) on the null simulations (R1-M4): CHM13 Hi-C (HiC-Pro, all MAPQ as in the
# paper) and HG002 Pore-C / CiFi (pairtools parse2 --min-mapq 1 as in the paper, and 0), at 50 kb and 20 kb.
set -euxo pipefail
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC
for r in 50000 20000; do
  bash $Q/contact_stats.sbatch --sample CHM13 --matrix $Q/chm13_null_pk40/sim.mcool --resolution $r --fai $HC/CHM13/CHM13.fa.gz.fai --out $Q/chm13_null_pk40/stats_$r
  for p in porec cifi; do for mq in 1 0; do
    bash $Q/contact_stats.sbatch --sample HG002 --matrix $Q/${p}_null_rk20/sim.mq$mq.mcool --resolution $r --fai $HC/HG002/HG002.fa.gz.fai --out $Q/${p}_null_rk20/stats_mq${mq}_$r
  done; done
done
