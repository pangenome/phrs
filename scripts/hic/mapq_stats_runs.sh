#!/bin/bash
# Run contact_stats.sbatch (as a plain script) on: the paper HG002 Hi-C 50 kb (reproduction with the rebuilt env), then
# the MAPQ >= 1 and >= 2 matrices of HG002 and CHM13 at 50 kb and 20 kb; then compare the reproduction with the paper.
set -euxo pipefail
Q=/moosefs/guarracino/HPRCv2/PHR_III/simhic; HC=/moosefs/guarracino/HPRCv2/PHR_III/HiC; A=/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based
bash $Q/contact_stats.sbatch --sample HG002 --matrix $HC/HG002/hicpro_output/cool/HG002.mcool --resolution 50000 --fai $HC/HG002/HG002.fa.gz.fai --out $Q/repro_hg002_hic_50kb_rebuilt
for f in hg002_global_test.tsv hg002_contact_matrix.tsv hg002_community_bootstrap_tests.tsv hg002_phr_pair_correlation.tsv; do cmp $A/50000bp/$f $Q/repro_hg002_hic_50kb_rebuilt/$f; echo "identical $f"; done
for g in HG002 CHM13; do for t in 1 2; do for r in 50000 20000; do
  bash $Q/contact_stats.sbatch --sample $g --matrix $Q/mapq/$g.mapq$t.mcool --resolution $r --fai $HC/$g/$g.fa.gz.fai --out $Q/mapq/stats_${g}_mapq${t}_${r}
done; done; done
